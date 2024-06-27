/*

Copyright (c) 2005-2024, University of Oxford.
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

#include <boost/assign.hpp>
#include <cxxtest/TestSuite.h>

#include "AbstractDrugDataStructure.hpp"
#include "BayesianInferer.hpp"
#include "LogisticDistribution.hpp"

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
        for (unsigned i = 0; i < samples.size(); i++)
        {
            temp_mean += samples[i];
            samples[i] = AbstractDataStructure::ConvertPic50ToIc50(samples[i]);
        }
        temp_mean /= samples.size();

        std::cout << "Original iC50 = " << AbstractDataStructure::ConvertPic50ToIc50(rData[0]) << ", mean of inferred samples = " << AbstractDataStructure::ConvertPic50ToIc50(temp_mean) << std::endl
                  << std::flush;

        return samples;
    }

public:
    void TestPIC50Inference()
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

        std::vector<double> samples_should_be = boost::assign::list_of(4.3655)(4.4255)(4.6057)(4.8525)(4.4393)(4.8862)(4.3602)(4.8596)(4.1970)(4.4684);

        // First we test the method that returns a whole vector of samples
        std::vector<double> samples = inferer_object.GetSampleMedianValue(10u);

        // Re-seed so we get the same answer getting them one at once
        RandomNumberGenerator::Instance()->Reseed(0);

        // Check these against the ones we should get.
        for (unsigned i = 0; i < samples_should_be.size(); i++)
        {
            double sample = inferer_object.GetSampleMedianValue();
            TS_ASSERT_DELTA(sample, samples[i], 1e-12);
            TS_ASSERT_DELTA(sample, samples_should_be[i], 1e-4);
        }

        // Coverage
        TS_ASSERT_THROWS_THIS(BayesianInferer inferer_object(TESTING),
                              "No known distribution for this parameter.");
    }

    void TestHillInference()
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

    void TestRepeatedCalls()
    {
        // To deal with a troublesome case in the TestTqtCompounds.hpp
        const double sigma_iks_pic50 = 0.139736283;
        const double duloxetine_iks_ic50 = 9.577; // pIc50 approx 5.2

        unsigned num_samples = 100;
        std::vector<double> samples;

        std::vector<double> data;
        data.push_back(AbstractDataStructure::ConvertIc50ToPic50(duloxetine_iks_ic50));

        for (unsigned i = 0; i < 10; i++)
        {
            samples = GetIc50Samples(data,
                                     sigma_iks_pic50,
                                     num_samples);
        }
        TS_ASSERT_EQUALS(samples.size(), num_samples);
    }

    // A long running test that triggered a refactor to use log-likelihoods. Don't run all the time.
    void xTestATroublesomeCase()
    {
        std::vector<double> large_pic50_dataset{ 5.82, 5.96, 6.03, 6.15, 5.99, 5.98, 6.03, 6.07, 6.16, 5.96, 6.03, 6.35, 6.72, 5.92, 6.05, 4.15, 5.9, 6.64, 6.13, 6.15, 6.07, 6.17, 6.54, 6.32, 6.23, 5.6, 5.88, 5.99, 6.24, 5.96, 6.03, 6.75, 5.92, 6.409, 5.53, 6.14, 6.37, 5.96, 5.6, 5.81, 6.23, 6.11, 5.54, 5.65, 6.4, 6.29, 6.37, 6.42, 6.09, 6.14, 6.05, 6.41, 5.77, 6.46, 5.91, 6.24, 5.63, 6.35, 5.86, 6.36, 6.15, 6.86, 6.11, 6.326, 6.16, 5.9, 6.03, 6.19, 6.44, 6.45, 6.14, 6.34, 6.7, 6.159, 6.14, 6.29, 5.6, 6.034, 6.02, 5.98, 5.92, 6.19, 5.55, 6.47, 6.11, 6.27, 6.39, 5.91, 4.3, 6.16, 6.21, 6.2, 5.74, 5.69, 6.08, 6.24, 6.69, 6.517, 5.57, 6.08, 6.01, 6.17, 6.29, 5.96, 6.25, 5.93, 5.87, 6.03, 6.42, 6.03, 5.95, 5.87, 6.92, 5.41, 6.01, 6.24, 6.71, 6.09, 6.25, 6.12, 6.21, 5.93, 5.99, 6.203, 6.107, 6.28, 6.22, 5.81, 6.21, 6.13, 6.2, 6, 6.46, 5.83, 6.27, 5.97, 5.55, 6.23, 6.25, 6.54, 6, 6.25, 6.21, 5.98, 5.86, 5.94, 5.98, 5.88, 5.98, 6.448, 5.79, 5.56, 6.3, 6.16, 6.19, 6.09, 6.11, 6.57, 6.108, 6.26, 5.99, 5.83, 5.86, 7.23, 6.2, 6.59, 5.79, 5.99, 6.04, 5.71, 6.21, 5.67, 6.03, 6.4, 6.79, 5.85, 6.03, 6.13, 5.99, 6.05, 5.84, 6.19, 6.17, 6.3, 6.427, 6.25, 5.94, 6.2, 6.29, 5.49, 6.86, 5.98, 6.39, 6.02, 6.02, 6.2, 6.03, 6.07, 6.18, 6.34, 6.51, 6.11, 6.1, 6.06, 6.52, 6.09, 5.6, 6.414, 6.77, 6.68, 6.12, 6.23, 6.1, 6.89, 6.13, 6.153, 5.98, 6.34, 6.1, 5.9, 5.62, 5.67, 5.96, 6.34, 6.019, 5.83, 6.317, 6.44, 5.93, 5.83, 5.96, 6.011, 6.52, 5.69, 6.04, 6.18, 6.25, 5.95, 6.33, 5.93, 5.88, 6.45, 5.47, 6.12, 6.08, 6.33, 6, 6.38, 6.58, 5.7, 6.81, 6.06, 6.34, 6.2, 6.28, 6.03, 5.91, 5.82, 6.36, 6.24, 6.046, 6.454, 6.37, 6.32, 5.98, 6.51, 6.35, 6.38, 6.25, 6.48, 6.13, 6.45, 6.18, 6.33, 6.23, 5.72, 6.133, 6.29, 6.19, 6.34, 6, 5.99, 5.93, 6.43, 6.48, 6.13, 6.369, 5.91, 6.135, 6.057, 6.02, 6.4, 5.99, 5.89, 6.2, 5.88, 5.71, 5.48, 6.41, 6.04, 6.1, 6.09, 6.09, 5.84, 5.97, 6.1, 5.72, 5.95, 6.29, 6.13, 6.3, 6.06, 5.95, 6.25, 6.25, 5.945, 5.86, 6.07, 5.9, 6.24, 5.88, 6.11, 6.4, 6.67, 6.19, 5.91, 6.032, 6.19, 5.97, 6.35, 6.15, 7.03, 5.94, 6.26, 6.09, 6.11, 6.09, 5.56, 6.09, 5.87, 5.91, 6.23, 6.4, 5.67, 6.18, 5.66, 5.86, 6.11, 6.44, 5.67, 5.99, 6.08, 6.12, 5.94, 6.54, 6.388, 5.73, 6.027, 6.02, 6.08, 6.08, 5.63, 6.05, 6.044, 6.07, 5.89, 6.28, 5.77, 6.34, 6.2, 6.29, 6.06, 5.49, 6.19, 6.36, 5.93, 6.15, 6.26, 6.32, 6.18, 6.22, 6, 6.37, 6.05, 6.4, 6.23, 5.93, 6.33, 6.26, 6.02, 6.06, 6.33, 6.2, 6.21, 6.26, 6.265, 5.98, 5.63, 6.29, 5.96, 6.49, 6.2, 6.23, 6.235, 5.91, 6.16, 6.12, 5.83, 6.77, 6.022, 6.01, 6.05, 6, 6.2, 5.94, 6.12, 6.14, 6.23, 6.45, 5.98, 6.19, 6.41, 5.96, 5.56, 5.8, 6.18, 6.4, 6.3, 5.77, 5.73, 6.03, 6.23, 6.12, 5.84, 5.88, 6.49, 6.04, 5.88, 6.33, 6.01, 6.14, 6.44, 5.99, 6.08, 6.08, 6.26, 6.29, 6.79, 5.97, 6.13, 5.84, 6.72, 6.56, 6.496, 5.87, 6.01, 5.6, 6.05, 6.13, 5.8, 6.39, 6.47, 6.37, 6.46, 6.27, 5.75, 6.09, 6.38, 6.1, 6.87, 6.31, 5.81, 6.16, 6.457, 5.963, 6.38, 6.33, 6.16, 6.31, 6.33, 6.03, 6.08, 6.17, 5.969, 6, 6.37, 5.93, 6.017, 5.79, 6.45, 6.21, 5.93, 6.65, 6.44, 6.06, 6.66, 6.27, 6.36, 6.43, 6.13, 5.94, 6.42, 6.28, 6.3, 6.3, 6.12, 6.27, 6.11, 5.92, 5.6, 6.341, 6.06, 6.3, 6.1, 6.15, 5.59, 5.62, 6.395, 5.88, 6.16, 6.25, 5.59, 6.166, 6.39, 6.112, 6.29, 6.28, 6.1, 5.83, 6.39, 6.35, 6.37, 6.2, 5.83, 6.162, 5.68, 6.53, 6.08, 5.88, 5.79, 6.03, 6.44, 5.91, 5.85, 6.05, 6.34, 6.31, 6, 5.97, 6.26, 6.07, 6.3, 6.1, 5.93, 6.14, 6.39, 6.036, 6.53, 5.76, 6.67, 6.35, 5.53, 5.918, 5.94, 6.16, 6.09, 6.06, 6.04, 6.29, 6.05, 5.81, 6.607, 5.99, 6.5, 5.79, 6.01, 6.11, 6.35, 6.12, 6.01, 5.9, 6.34, 5.57, 6.1, 6.42, 6.55, 6.71, 6.23, 5.94, 5.98, 6.2, 5.62, 6.3, 6.21, 6.42, 5.75, 6.18, 6.18, 6.15, 6.093, 6.25, 6.41, 6.71, 6.18, 5.48, 6.497, 5.8, 5.92, 5.87, 6.1, 6.44, 5.61, 6.14, 6.31, 6.07, 6.281, 6.15, 5.88, 6.14, 6.17, 6.53, 6.08, 6.66, 6.28, 5.98, 6.21, 6.63, 5.76, 6.06, 5.93, 6.36, 6.17, 5.91, 6.26, 5.43, 6.18, 6.35, 5.62, 6.19, 5.71, 5.54, 5.71, 6.54, 6.59, 5.9, 5.964, 6.8, 5.84, 5.82, 6.35, 6.11, 6.26, 5.66, 6.7, 6.47, 6.41, 6.18, 6.45, 5.9, 5.79, 6.64, 6, 6.27, 6.26, 6.54, 6.23, 6.13, 7.02, 6.155, 6.51, 6.45, 6.47, 6.02, 6.507, 5.84, 5.96, 6.19, 6.33, 6.36, 6.22, 6.391, 5.61, 6.05, 6.43, 5.89, 5.7, 5.89, 6.05, 5.89, 6.46, 5.77, 5.99, 6.18, 6.39, 6.33, 5.93, 5.64, 5.76, 6.36, 6.3, 5.92, 6.29, 6.08, 5.4013923292191119429617174318991601467132568359375, 5.63773038924257985371468748780898749828338623046875, 5.48368506831286328662145024281926453113555908203125, 4.7433987307798179955398154561407864093780517578125, 0, 0, 0, 4.583254542524077379539448884315788745880126953125, 5.35506565248765920017604003078304231166839599609375, 5.10899483732673420632863781065680086612701416015625, 4.848682378566724793245157343335449695587158203125, 0, 6.17293134723838665678385950741358101367950439453125, 5.42952944856429109421469547669403254985809326171875, 6.12105043985217367463746995781548321247100830078125, 0, 0, 4.3223052210139645268327512894757091999053955078125, 0, 4.633027433782682891205695341341197490692138671875, 4.76685028827561296793646761216223239898681640625, 0, 0, 0, 4.7144456163016972283230643370188772678375244140625, 4.80306298746282767098136901040561497211456298828125, 0, 0, 4.973252868491808698081513284705579280853271484375, 0, 5.05156864556507567698417915380559861660003662109375, 4.91489902537812728411381613113917410373687744140625, 0, 4.67111258940422491292565609910525381565093994140625, 0, 0, 5.30640274338680750787489159847609698772430419921875, 5.3446481215553571786358588724397122859954833984375, 5.3621623738753640964205260388553142547607421875, 5.219190190624264147345456876792013645172119140625, 5.45466120965008460785838906303979456424713134765625, 5.804184171748591580808351864106953144073486328125, 5.98677166477913313968883812776766717433929443359375, 5.7146542641841762133481097407639026641845703125, 0, 7.42700068290905246470856582163833081722259521484375, 5.3446481215553571786358588724397122859954833984375, 0, 4.999839167108316217991159646771848201751708984375, 5.804184171748591580808351864106953144073486328125, 8.4009907601766560247824600082822144031524658203125, 6.16731740423529917194400695734657347202301025390625, 5.696311373403258215830646804533898830413818359375, 6.472756883611911149500883766449987888336181640625, 5.54094271182182129820148475118912756443023681640625, 5.8915455933872973304232800728641450405120849609375, 5.1629529365534523321912274695932865142822265625, 5.2791880909903863283716418663971126079559326171875, 5.8392836285950266983491019345819950103759765625, 7.98653076769082925778775461367331445217132568359375, 4.528708288941061255172826349735260009765625, 5.3621623738753640964205260388553142547607421875, 5.0833436384358439141806229599751532077789306640625, 5.908919648029570481639893841929733753204345703125, 5, 9, 5, 6, 5, 5, 6, 6, 6, 5, 6, 6, 5, 5, 6, 6, 5, 6, 5, 6, 6, 8, 5, 6, 6, 5, 6, 6, 5, 6, 5, 6, 6, 6, 5, 6, 6, 6, 5, 6, 6, 5, 6, 6, 5, 5, 5, 5, 5, 5, 6, 5, 5, 5, 6, 5, 6, 6, 5 };
        BayesianInferer inferer_object(PIC50);

        double sigma = 0.5;

        inferer_object.SetObservedData(large_pic50_dataset);
        inferer_object.SetSpreadOfUnderlyingDistribution(sigma);

        inferer_object.PerformInference();

        std::vector<double> values = inferer_object.GetPossibleMedianValues();
        std::vector<double> pdf = inferer_object.GetPosteriorPdf();

        for (unsigned i = 0; i < values.size(); i++)
        {
            std::cout << values[i] << "\t" << pdf[i] << std::endl;
        }

        TS_ASSERT_DELTA(inferer_object.GetSampleMedianValue(), 6.0059, 1e-3);
    }
};

#endif // TESTBAYESIANINFERER_HPP_
