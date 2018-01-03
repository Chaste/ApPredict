/*

Copyright (c) 2005-2018, University of Oxford.
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

#ifndef TESTBAYESIANINFERER_HPP_
#define TESTBAYESIANINFERER_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/assign.hpp>

#include "BayesianInferer.hpp"
#include "LogisticDistribution.hpp"
#include "AbstractDrugDataStructure.hpp"

class TestBayesianInferer : public CxxTest::TestSuite
{
private:
    std::vector<double> GetIc50Samples(std::vector<double>& rData,
                                       const double& rSigma,
                                       const unsigned& rNumSamples)
    {
        // We will say that the underlying distribution really is
        // a logistic with this ic50 and sigma.
        LogisticDistribution logistic;
        BayesianInferer inferer(PIC50);

        inferer.SetObservedData(rData);
        inferer.SetSpreadOfUnderlyingDistribution(rSigma);
        inferer.PerformInference(); // Generate a PDF and CDF for the underlying distribution centring parameter.

        // We now take samples from the inferred PDF
        std::vector<double> samples = inferer.GetSampleMedianValue(rNumSamples);

        // Convert them from pIC50 to IC50 for dose-response calculations.
        double temp_mean = 0u;
        for (unsigned i=0; i<samples.size(); i++)
        {
            temp_mean += samples[i];
            samples[i] = AbstractDataStructure::ConvertPic50ToIc50(samples[i]);
        }
        temp_mean /= samples.size();

        std::cout << "Original iC50 = " << AbstractDataStructure::ConvertPic50ToIc50(rData[0]) <<
                ", mean of inferred samples = " << AbstractDataStructure::ConvertPic50ToIc50(temp_mean) << std::endl << std::flush;

        return samples;
    }

public:
    void TestPIC50Inference() throw (Exception)
    {
        BayesianInferer inferer_object(PIC50);

        TS_ASSERT_THROWS_THIS(inferer_object.GetSampleMedianValue(),
                "Inference has not been performed, please call PerformInference() before trying to get samples.");

        TS_ASSERT_THROWS_THIS(inferer_object.PerformInference(),
                "Please call SetObservedData() and SetSpreadOfUnderlyingDistribution() before PerformInference().");

        std::vector<double> data;
        data.push_back(4.2);
        data.push_back(4.4);

        double sigma = 0.5;

        inferer_object.SetObservedData(data);
        inferer_object.SetSpreadOfUnderlyingDistribution(sigma);

        inferer_object.PerformInference();

        std::vector<double> samples_should_be = boost::assign::list_of(4.3655)(4.4255)(4.6057)(4.8525)
                (4.4393)(4.8862)(4.3602)(4.8596)(4.1970)(4.4684);

        // First we test the method that returns a whole vector of samples
        std::vector<double> samples = inferer_object.GetSampleMedianValue(10u);

        // Re-seed so we get the same answer getting them one at once
        RandomNumberGenerator::Instance()->Reseed(0);

        // Check these against the ones we should get.
        for (unsigned i=0; i<samples_should_be.size(); i++)
        {
            double sample = inferer_object.GetSampleMedianValue();
            TS_ASSERT_DELTA(sample, samples[i], 1e-12);
            TS_ASSERT_DELTA(sample, samples_should_be[i], 1e-4);
        }

        // Coverage
        TS_ASSERT_THROWS_THIS(BayesianInferer inferer_object(TESTING),
                "No known distribution for this parameter.");
    }

    void TestHillInference() throw (Exception)
    {
        BayesianInferer inferer_object(HILL);

        std::vector<double> data;
        data.push_back(1.2);
        data.push_back(0.9);

        double beta = 4.1;

        inferer_object.SetObservedData(data);
        inferer_object.SetSpreadOfUnderlyingDistribution(beta);

        TS_ASSERT_THROWS_THIS(inferer_object.GetPosteriorCdf(),
                "Posterior has not yet been computed, call PerformInference() first.");

        inferer_object.PerformInference();

        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 1.2414, 1e-4);
        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 1.0351, 1e-4);
        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 1.0739, 1e-4);
        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 0.9711, 1e-4);
        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 1.5963, 1e-4);
        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 0.7255, 1e-4);
        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 1.9176, 1e-4);
        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 0.9521, 1e-4);
        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 1.0344, 1e-4);
        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 1.1036, 1e-4);

        std::vector<double> posterior = inferer_object.GetPosteriorCdf();
        std::vector<double> values = inferer_object.GetPossibleMedianValues();
        TS_ASSERT_EQUALS(posterior.size(), values.size());

        // Following snippet is useful for copying and pasting output into matlab for plotting.
//        for(unsigned i=0; i<values.size(); i++)
//        {
//            std::cout << values[i] << "\t" << posterior[i] << "\n";
//        }
    }

    void TestRepeatedCalls() throw (Exception)
    {
        // To deal with a troublesome case in the TestTqtCompounds.hpp
        const double sigma_iks_pic50 = 0.139736283;
        const double duloxetine_iks_ic50 = 9.577; // pIc50 approx 5.2

        unsigned num_samples = 100;
        std::vector<double> samples;

        std::vector<double> data;
        data.push_back(AbstractDataStructure::ConvertIc50ToPic50(duloxetine_iks_ic50));

        for (unsigned i=0; i<10; i++)
        {
            samples = GetIc50Samples(data,
                                     sigma_iks_pic50,
                                     num_samples);
        }
        TS_ASSERT_EQUALS(samples.size(),num_samples);
    }
};

#endif // TESTBAYESIANINFERER_HPP_
