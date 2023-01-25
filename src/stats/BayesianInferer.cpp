/*

Copyright (c) 2005-2023, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <cmath>
#include <vector>

#include "Exception.hpp"

#include "BayesianInferer.hpp"
#include "LogisticDistribution.hpp"
#include "LogLogisticDistribution.hpp"

// Add the citation for Elkins paper
#include "Citations.hpp"
static PetscBool ElkinsCite = PETSC_FALSE;
const char ElkinsCitation[] = "@article{elkins2013variability,\n"
                              "  title={Variability in high-throughput ion-channel screening data and "
                              "consequences for cardiac safety assessment},\n"
                              "  author={Elkins, Ryan C and Davies, Mark R and Brough, Stephen J and "
                              "Gavaghan, David J and Cui, Yi and Abi-Gerges, Najah and Mirams, Gary R},\n"
                              "  journal={Journal of pharmacological and toxicological methods},\n"
                              "  volume={68},\n"
                              " number={1},\n"
                              "  pages={112--122},\n"
                              "  year={2013},\n"
                              "  publisher={Elsevier}\n"
                              "}";

BayesianInferer::BayesianInferer(DoseResponseParameter parameter)
    : mParameter(parameter),
      mSigma(DOUBLE_UNSET),
      mInferenceReady(false),
      mpData(NULL)
{
    // Record a reference for the calculations performed here, can be extracted with the '-citations' flag.
    Citations::Register(ElkinsCitation, &ElkinsCite);

    unsigned num_values = 1000000u; // it is very quick so put loads and loads of points in for nice smooth distributions.
    double max_value;
    double min_value;

    // Set up the limits on possible mu values for this parameter (to form boundaries of our prior).
    if (mParameter == PIC50)
    {
        max_value = 12;
        min_value = -12;
        mpDistribution = new LogisticDistribution();
    }
    else if (mParameter == HILL)
    {
        min_value = 0.1;
        max_value = 10;
        mpDistribution = new LogLogisticDistribution();
    }
    else
    {
        EXCEPTION("No known distribution for this parameter.");
    }

    // Set up a range of possible mu values for this parameter
    for (unsigned i = 0; i < num_values; i++)
    {
        mPossibleMuValues.push_back(min_value + ((double)(i)) * (max_value - min_value) / ((double)(num_values)-1.0));
    }
};

BayesianInferer::~BayesianInferer()
{
    if (mpDistribution)
    {
        delete mpDistribution;
    }
};

void BayesianInferer::SetObservedData(const std::vector<double> &rData)
{
    // Point directly to the 'const' data.
    mpData = &rData;
}

void BayesianInferer::SetSpreadOfUnderlyingDistribution(double sigma)
{
    mSigma = sigma;
}

void BayesianInferer::PerformInference()
{
    if (mSigma == DOUBLE_UNSET || mpData == NULL)
    {
        EXCEPTION("Please call SetObservedData() and SetSpreadOfUnderlyingDistribution() before PerformInference().");
    }

    // Set up our prior distribution - it is uniform, so just divide one by number of possible options.
    unsigned num_possible_values = mPossibleMuValues.size();
    const double log_prior_prob_this_mu = log(1.0 / ((double)(num_possible_values)));
    mPosteriorPdf.resize(num_possible_values);
    mPosteriorCdf.resize(num_possible_values);
    std::vector<double> log_posterior(num_possible_values);

    // Calculate the modifications to the prior distribution
    for (unsigned i = 0; i < num_possible_values; i++)
    {
        log_posterior[i] = log_prior_prob_this_mu;

        for (unsigned j = 0; j < mpData->size(); j++)
        {
            log_posterior[i] += log(mpDistribution->EvaluatePdf(mPossibleMuValues[i], mSigma, (*mpData)[j]));
        }
    }

    double max_log_likelihood = *std::max_element(log_posterior.cbegin(), log_posterior.cend());

    // Do some scaling so that our posterior distribution is a PDF.
    // First, work out the sum of the posterior likelihoods (scaled up to sensible numbers).
    for (unsigned i = 0; i < num_possible_values; i++)
    {
        mPosteriorPdf[i] = exp(log_posterior[i] - max_log_likelihood); // Add a constant to get these order one and nice for computations.
    }
    double sum = std::accumulate(mPosteriorPdf.begin(), mPosteriorPdf.end(), 0.0);
    assert(sum > 0.0);

    double scaling_factor = (((double)(num_possible_values)) / (mPossibleMuValues.back() - mPossibleMuValues[0])) / sum;
    for (unsigned i = 0; i < num_possible_values; i++)
    {
        mPosteriorPdf[i] *= scaling_factor;
    }

    // Do a simple sum to work out the CDF from the PDF
    mPosteriorCdf[0] = mPosteriorPdf[0];
    for (unsigned i = 1; i < num_possible_values; i++)
    {
        mPosteriorCdf[i] = mPosteriorCdf[i - 1u] + mPosteriorPdf[i];
    }
    // And scale so that it really is a CDF.
    for (unsigned i = 0; i < num_possible_values; i++)
    {
        mPosteriorCdf[i] /= mPosteriorCdf.back();
    }

    // Safety checks.
    assert(mPosteriorCdf.size() == mPossibleMuValues.size());
    assert(mPosteriorPdf.size() == mPossibleMuValues.size());

    mInferenceReady = true;
}

double BayesianInferer::GetSampleMedianValue()
{
    if (!mInferenceReady)
    {
        EXCEPTION("Inference has not been performed, please call PerformInference() before trying to get samples.");
    }
    assert(mPosteriorCdf.size() == mPossibleMuValues.size());

    // Get a random number p in [0,1] to use as backwards lookup in posterior CDF.
    double p = RandomNumberGenerator::Instance()->ranf();

    // Loop through the CDF to see when it goes above this random number `p',
    // then do a bit of interpolation to give a unique sample 'x' for this `p'.
    unsigned index = UNSIGNED_UNSET;
    for (unsigned i = 0; i < mPosteriorCdf.size(); i++)
    {
        if (mPosteriorCdf[i] >= p)
        {
            index = i;
            break;
        }
    }
    assert(index != UNSIGNED_UNSET);
    double return_value = mPossibleMuValues[index];
    if (index < mPosteriorCdf.size() - 1u)
    {
        // Do a linear interpolation to ensure unique answers each call.
        double proportion_through = (p - mPosteriorCdf[index]) / (mPosteriorCdf[index + 1] - mPosteriorCdf[index]);
        return_value += proportion_through * (mPossibleMuValues[index + 1] - mPossibleMuValues[index]);
    }
    return return_value;
}

std::vector<double> BayesianInferer::GetSampleMedianValue(const unsigned numValues)
{
    std::vector<double> samples(numValues);
    for (unsigned i = 0; i < numValues; i++)
    {
        samples[i] = GetSampleMedianValue();
    }

    assert(samples.size() == numValues);
    return samples;
}

std::vector<double> BayesianInferer::GetPossibleMedianValues()
{
    return mPossibleMuValues;
}

std::vector<double> BayesianInferer::GetPosteriorCdf()
{
    if (!mInferenceReady)
    {
        EXCEPTION("Posterior has not yet been computed, call PerformInference() first.");
    }
    return mPosteriorCdf;
}

std::vector<double> BayesianInferer::GetPosteriorPdf()
{
    if (!mInferenceReady)
    {
        EXCEPTION("Posterior has not yet been computed, call PerformInference() first.");
    }
    return mPosteriorPdf;
}

double BayesianInferer::GetSpreadOfUnderlyingDistribution()
{
    if (!mInferenceReady)
    {
        EXCEPTION("Posterior has not yet been computed, call PerformInference() first.");
    }
    return mSigma;
}