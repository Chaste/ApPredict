/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTLOGISTICDISTRIBUTION_HPP_
#define TESTLOGISTICDISTRIBUTION_HPP_

#include <cxxtest/TestSuite.h>
#include <numeric>
#include <vector>

#include "LogLogisticDistribution.hpp"
#include "LogisticDistribution.hpp"

class TestLogisticDistribution : public CxxTest::TestSuite
{
private:
    static const unsigned mNumRepeats = 1000000;

public:
    void TestBasicDistribution()
    {
        double mu = 2.0;
        double sigma = 1.0;
        LogisticDistribution sampler;

        std::vector<double> values;
        for (unsigned i = 0; i < mNumRepeats; i++)
        {
            values.push_back(sampler.GetSample(mu, sigma));
        } // End of loop over samples

        // Check that the mean and standard deviation of the values is correct
        double theoretical_mean = mu;
        double theoretical_std = sigma * M_PI / sqrt(3.0);

        // Calculate mean and std of 'values'
        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        double recorded_mean = sum / values.size();
        double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
        double recorded_std = std::sqrt(sq_sum / values.size() - recorded_mean * recorded_mean);

        TS_ASSERT_DELTA(theoretical_mean, recorded_mean, 2e-3);
        TS_ASSERT_DELTA(theoretical_std, recorded_std, 2e-3);
    }

    /**
     * Here we test the probability density function (PDF).
     */
    void TestPdfCalculations()
    {
        double mu = 4.5;
        double sigma = 0.3;
        double sample = 5.0;
        LogisticDistribution distribution;
        double probability = distribution.EvaluatePdf(mu, sigma, sample);

        TS_ASSERT_DELTA(probability, 0.44543237465084, 1e-12);

        // There is no really nice relationship to compare easily
        // as the PDF depends on the x-scale which is very different.
        LogLogisticDistribution distribution2;
        probability = distribution2.EvaluatePdf(1, 8, 1);

        TS_ASSERT_DELTA(probability, 2.0, 1e-12);
    }

    /**
     * Samples are calculated using the inverse CDF or 'Quantile' function.
     */
    void TestMultipleSamples()
    {
        double mu = 2.0;
        double sigma = 1.0;
        unsigned num_experiments = 4;
        LogisticDistribution sampler;

        std::vector<double> values;

        for (unsigned i = 0; i < mNumRepeats; i++)
        {
            values.push_back(sampler.GetSample(mu, sigma, num_experiments));
        } // End of loop over samples

        // Check that the mean and standard deviation of the values is correct
        double theoretical_mean = mu;
        double theoretical_std = sigma * M_PI / sqrt(3.0 * num_experiments);

        // Calculate mean and std of 'values'
        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        double recorded_mean = sum / values.size();
        double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
        double recorded_std = std::sqrt(sq_sum / values.size() - recorded_mean * recorded_mean);

        TS_ASSERT_DELTA(theoretical_mean, recorded_mean, 2e-3);
        TS_ASSERT_DELTA(theoretical_std, recorded_std, 2e-3);
    }

    /**
     * Samples are calculated using the inverse CDF or 'Quantile' function.
     */
    void TestLogLogisticSamples()
    {
        double alpha = exp(-0.5);
        double beta = 1.0 / 0.12;
        unsigned num_experiments = 1;
        LogLogisticDistribution sampler;

        std::vector<double> values;

        for (unsigned i = 0; i < mNumRepeats; i++)
        {
            values.push_back(sampler.GetSample(alpha, beta, num_experiments));
        } // End of loop over samples

        // Check that the mean and standard deviation of the values is correct
        double b = M_PI / beta;
        double theoretical_mean = alpha * b / sin(b);
        double theoretical_variance = alpha * alpha * ((2 * b) / (sin(2 * b)) - b * b / (sin(b) * sin(b)));

        // Calculate mean and std of 'values'
        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        double recorded_mean = sum / values.size();
        double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
        double recorded_variance = sq_sum / values.size() - recorded_mean * recorded_mean;

        TS_ASSERT_DELTA(theoretical_mean, recorded_mean, 1e-3);
        TS_ASSERT_DELTA(theoretical_variance, recorded_variance, 1e-3);
    }
};

#endif // TESTLOGISTICDISTRIBUTION_HPP_
