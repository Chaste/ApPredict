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
#ifdef CHASTE_CVODE

#ifndef _TESTAPPREDICTLONG_HPP_
#define _TESTAPPREDICTLONG_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/shared_ptr.hpp>
#include "ApPredictMethods.hpp"
#include "CommandLineArgumentsMocker.hpp"
#include "NumericFileComparison.hpp"

/*
 * Thorough and longer running tests of ApPredict.
 *
 * Note we redirect output from ApPredict to avoid output folder conflict if running this test in parallel.
 */
class TestApPredictLong : public CxxTest::TestSuite
{
public:
    /**
     * This test will wipe $CHASTE_TEST_OUTPUT/ApPredict_output_long/
     *
     * The tests overwrite CommandLineArguments and does some standard
     * simulations to check things are working OK...
     */
    void TestSomeFavouriteCompounds(void)
    {
        // Test a simple hERG block with TT06
        // loop over hardcoded and dynamically loaded.
        // try both --model --cellml and --cellml and try both relative and absolute paths

        FileFinder cellm_file("projects/ApPredict/src/cellml/ten_tusscher_model_2006_epi.cellml", RelativeTo::CWD);
        std::vector<std::string> model_args = { "--model 2",
                                                "--model ten_tusscher_model_2006_epi",
                                                "--model projects/ApPredict/src/cellml/ten_tusscher_model_2006_epi.cellml",
                                                "--model " + cellm_file.GetAbsolutePath(),
                                                "--cellml projects/ApPredict/src/cellml/ten_tusscher_model_2006_epi.cellml"};
        for (unsigned i = 0; i < model_args.size(); i++)
        {
            CommandLineArgumentsMocker wrapper(model_args[i] + " --plasma-concs 1 10 --pic50-herg 4.5 --plasma-conc-logscale false --output-dir ApPredict_output_long");

            ApPredictMethods methods;
            methods.Run();

            std::vector<double> concs = methods.GetConcentrations();

            TS_ASSERT_EQUALS(concs.size(), 3u);
            TS_ASSERT_DELTA(concs[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(concs[1], 1.0, 1e-12);
            TS_ASSERT_DELTA(concs[2], 10.0, 1e-12);

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(), 3u);
            TS_ASSERT_DELTA(apd90s[0], 301.4617, 1e-3);
            TS_ASSERT_DELTA(apd90s[1], 302.6703, 1e-3);
            TS_ASSERT_DELTA(apd90s[2], 311.27, 1e-2);

            TS_ASSERT_THROWS_THIS(methods.GetApd90CredibleRegions(),
                                  "There was no Lookup Table available for credible interval calculations with these settings.");
        }
    }

    void TestWithConfidenceIntervalOutputs(void)
    {
        // Test a couple of Exceptions
        {
            CommandLineArgumentsMocker wrapper("--model 1 --plasma-concs 1 10 --pic50-ik1 4.5 --credible-intervals --plasma-conc-logscale false --output-dir ApPredict_output_long");
            ApPredictMethods methods;
            TS_ASSERT_THROWS_THIS(methods.Run(),
                                  "No argument --pic50-spread-ik1 has been provided. Cannot calculate credible intervals without this.");
        }

        {
            CommandLineArgumentsMocker wrapper("--model 1 --plasma-concs 1 10 --pic50-herg 4.5 --credible-intervals --plasma-conc-logscale false --output-dir ApPredict_output_long");
            ApPredictMethods methods;
            TS_ASSERT_THROWS_THIS(methods.Run(),
                                  "No argument --pic50-spread-herg has been provided. Cannot calculate credible intervals without this.");
        }

        // Test a simple hERG block with Shannon
        {
            CommandLineArgumentsMocker wrapper("--model 1 --plasma-concs 1 10 --pic50-herg 4.5 --pic50-spread-herg 0.15 --credible-intervals --plasma-conc-logscale false --output-dir ApPredict_output_long");

            ApPredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();

            TS_ASSERT_EQUALS(concs.size(), 3u);
            TS_ASSERT_DELTA(concs[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(concs[1], 1.0, 1e-12);
            TS_ASSERT_DELTA(concs[2], 10.0, 1e-12);

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(), 3u);
            TS_ASSERT_DELTA(apd90s[0], 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90s[1], 213.9762, 1e-3);
            TS_ASSERT_DELTA(apd90s[2], 229.7436, 1e-3);

            std::vector<std::vector<double>> apd90_credible_regions = methods.GetApd90CredibleRegions();
            TS_ASSERT_EQUALS(apd90_credible_regions.size(), 3u);
            TS_ASSERT_DELTA(apd90_credible_regions[0][0], 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[0][1], 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[1][0], 212.5283, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[1][1], 219.0559, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[2][0], 217.5297, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[2][1], 259.4240, 1e-3);
        }

        // We should get much reduced credible regions with repeated pIC50 values
        {
            CommandLineArgumentsMocker wrapper("--model 1 --plasma-concs 1 10 --pic50-herg 4.5 4.4 4.6 4.5 --pic50-spread-herg 0.15 --credible-intervals --plasma-conc-logscale false --output-dir ApPredict_output_long");

            ApPredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();

            TS_ASSERT_EQUALS(concs.size(), 3u);
            TS_ASSERT_DELTA(concs[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(concs[1], 1.0, 1e-12);
            TS_ASSERT_DELTA(concs[2], 10.0, 1e-12);

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(), 3u);
            TS_ASSERT_DELTA(apd90s[0], 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90s[1], 213.9730, 1e-3);
            TS_ASSERT_DELTA(apd90s[2], 229.7191, 1e-3);

            std::vector<std::vector<double>> apd90_credible_regions = methods.GetApd90CredibleRegions();
            TS_ASSERT_EQUALS(apd90_credible_regions.size(), 3u);
            TS_ASSERT_DELTA(apd90_credible_regions[0][0], 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[0][1], 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[1][0], 213.1849, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[1][1], 215.3668, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[2][0], 223.2300, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[2][1], 239.4459, 1e-3);
        }

        // This is a case Gef found that seems to get 'mean' predictions outside the credible region
        // It turned out this was an error we had introduced by taking the mean IC50 value (instead of median PIC50),
        // or even better the median of the inferred pIC50s from the credible interval pIC50 samples.
        // (which is what we are doing now).
        {
            CommandLineArgumentsMocker wrapper("--model 1 --pacing-freq 0.5 "
                                               "--pic50-herg 4.86 4.45 0 4.1429 --hill-herg 3.22 1.57 1 1 "
                                               "--pic50-spread-herg 0.15 --hill-spread-herg 0.21 "
                                               "--credible-intervals --plasma-concs 30.0 --output-dir ApPredict_output_long");

            ApPredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();
            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(concs.size(), apd90s.size());
            TS_ASSERT_DELTA(apd90s[0], 219.604, 1e-1);
            TS_ASSERT_DELTA(apd90s[1], 219.604, 1e-1);
            TS_ASSERT_DELTA(apd90s[2], 247.4263, 1e-1);
        }

        // This is the same case as above but we say we have repeated IC50s but no spread info
        // In this case all we can do is take the median of the PIC50s. This gives a slightly
        // different hERG IC50 of 50.5301 uM
        {
            CommandLineArgumentsMocker wrapper("--model 1 --pacing-freq 0.5 "
                                               "--pic50-herg 4.86 4.45 0 4.1429 --hill-herg 3.22 1.57 1 1 "
                                               "--plasma-concs 30.0  --output-dir ApPredict_output_long");

            ApPredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();
            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(concs.size(), apd90s.size());
            TS_ASSERT_DELTA(apd90s[0], 219.604, 1e-1);
            TS_ASSERT_DELTA(apd90s[1], 219.604, 1e-1);
            TS_ASSERT_DELTA(apd90s[2], 250.522, 1e-1);
        }

        // This is a new case where we have an agonist (activator) rather than an
        // antagonist (inhibitor) of a channel.
        {
            CommandLineArgumentsMocker wrapper("--model 1 --pacing-freq 0.5 "
                                               "--pic50-herg 4.86 4.45 0 4.1429 --hill-herg 3.22 1.57 1 1 "
                                               "--plasma-concs 30.0 --saturation-herg 110 120 130.2 105.7 --output-dir ApPredict_output_long");

            ApPredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();
            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(concs.size(), apd90s.size());
            TS_ASSERT_DELTA(apd90s[0], 219.604, 1e-1);
            TS_ASSERT_DELTA(apd90s[1], 219.604, 1e-1);
            TS_ASSERT_DELTA(apd90s[2], 215.921, 1e-1); // hERG activator shortens AP.
        }
    }

    void TestWithAModelInAlternans(void)
    {
        // We are trying to go so fast we get alternans, and then handle it nicely.
        CommandLineArgumentsMocker wrapper("--model 3 --plasma-concs 1 10 --pacing-freq 5 --plasma-conc-logscale false --output-dir ApPredict_output_long");

        ApPredictMethods methods;
        methods.Run();
        std::vector<double> concs = methods.GetConcentrations();

        TS_ASSERT_EQUALS(concs.size(), 3u);
        TS_ASSERT_DELTA(concs[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(concs[1], 1.0, 1e-12);
        TS_ASSERT_DELTA(concs[2], 10.0, 1e-12);

        std::vector<double> apd90s = methods.GetApd90s();
        TS_ASSERT_EQUALS(apd90s.size(), 3u);
        TS_ASSERT_DELTA(apd90s[0], 137.3714, 1e-2);
        TS_ASSERT_DELTA(apd90s[1], 137.3714, 1e-2);
        TS_ASSERT_DELTA(apd90s[2], 137.3714, 1e-2);
    }

    void TestWithAModelTwoToOneStimAp(void)
    {
        // We should go 'too fast' for O'Hara and see if we can get anything like sensible output.
        CommandLineArgumentsMocker wrapper("--model 6 --plasma-concs 1 --pacing-freq 5 --plasma-conc-logscale false --output-dir ApPredict_output_long");

        ApPredictMethods methods;
        methods.Run();
        std::vector<double> concs = methods.GetConcentrations();

        TS_ASSERT_EQUALS(concs.size(), 2u);
        TS_ASSERT_DELTA(concs[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(concs[1], 1.0, 1e-12);

        std::vector<double> apd90s = methods.GetApd90s();
        TS_ASSERT_EQUALS(apd90s.size(), 2u);
        TS_ASSERT_DELTA(apd90s[0], 266.81, 1e-1); // Very sensitive to compiler changes, so high tolerance.
        TS_ASSERT_DELTA(apd90s[1], 266.81, 1e-1); // Very sensitive to compiler changes, so high tolerance.
    }

    void TestTroublesomeApCalculation(void)
    {
        // This is a case that Gef was having some trouble with. It does look as if something odd is going
        // on, as it reports AP alternans, then tries to swap round the pacing and it has gone, at 10k paces.
        CommandLineArgumentsMocker wrapper("--model 1 --pacing-freq 0.5 --pic50-cal 3.0 --hill-cal 1 "
                                           "--pic50-herg 0 --hill-herg 1 "
                                           "--pic50-na 4.561 --hill-na 1 "
                                           "--plasma-concs 0 100.0 --output-dir ApPredict_output_long");

        ApPredictMethods methods;
        methods.Run();
        std::vector<double> concs = methods.GetConcentrations();

        TS_ASSERT_EQUALS(concs.size(), 3u);

        std::vector<double> apd90s = methods.GetApd90s();
        TS_ASSERT_EQUALS(apd90s.size(), 3u);
        TS_ASSERT_DELTA(apd90s[0], 219.60, 1e-1); // Very sensitive to compiler changes, so high tolerance.
        TS_ASSERT_DELTA(apd90s[1], 219.60, 1e-1);
        TS_ASSERT_EQUALS(std::isnan(apd90s[2]), true); // Top concentration has a NoAP1 error.
    }

    void TestTwoDrugs()
    {
        // Check single drug version
        std::cout << "\nSingle drug:\n"
                  << std::endl;
        {
            CommandLineArgumentsMocker wrapper("--model 2 --plasma-concs 0 10 --pic50-herg 7 --pacing-max-time 0.2 --credible-intervals --pic50-spread-herg 0.15");

            ApPredictMethods methods;
            methods.SetOutputDirectory("ApPredict_only_one_drug/");
            methods.Run();
        }
        std::cout << "\nFirst drug potent:\n"
                  << std::endl;
        // Check for symmetry in the responses
        {
            CommandLineArgumentsMocker wrapper("--model 2 --plasma-concs 0 10 --pic50-herg 7 --drug-two-conc-factor 1 --pic50-drug-two-herg -10 --pacing-max-time 0.2 --credible-intervals --pic50-spread-herg 0.15 --pic50-spread-drug-two-herg 0.15");

            ApPredictMethods methods;
            methods.SetOutputDirectory("ApPredict_drug_1_potent/");
            methods.Run();
        }
        std::cout << "\nSecond drug potent:\n"
                  << std::endl;
        {
            CommandLineArgumentsMocker wrapper("--model 2 --plasma-concs 0 10 --pic50-herg -10 --drug-two-conc-factor 1 --pic50-drug-two-herg 7 --pacing-max-time 0.2 --credible-intervals --pic50-spread-herg 0.15 --pic50-spread-drug-two-herg 0.15");

            ApPredictMethods methods;
            methods.SetOutputDirectory("ApPredict_drug_2_potent/");
            methods.Run();
        }

        FileFinder drug_1_active_file("ApPredict_drug_1_potent/voltage_results.dat", RelativeTo::ChasteTestOutput);
        FileFinder drug_2_active_file("ApPredict_drug_2_potent/voltage_results.dat", RelativeTo::ChasteTestOutput);
        TS_ASSERT(drug_1_active_file.IsFile());
        TS_ASSERT(drug_2_active_file.IsFile());

        NumericFileComparison comparer(drug_1_active_file, drug_2_active_file);
        TS_ASSERT(comparer.CompareFiles(3e-2));
    }
};

#endif //_TESTAPPREDICTLONG_HPP_

#endif //_CHASTE_CVODE
