/*

Copyright (c) 2005-2015, University of Oxford.
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
#include "CommandLineArgumentsMocker.hpp"
#include "ApPredictMethods.hpp"

class TestApPredictLong : public CxxTest::TestSuite
{
public:

    /**
     * This test will wipe $CHASTE_TEST_OUTPUT/ApPredict_output/
     *
     * The tests overwrite CommandLineArguments and does some standard
     * simulations to check things are working OK...
     */
    void TestSomeFavouriteCompounds(void) throw (Exception)
    {
        // Test a simple hERG block with TT06
        {
            CommandLineArgumentsMocker wrapper("--model 2 --plasma-concs 1 10 --pic50-herg 4.5 --plasma-conc-logscale false");

            ApPredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();

            TS_ASSERT_EQUALS(concs.size(),3u);
            TS_ASSERT_DELTA(concs[0], 0.0,   1e-12);
            TS_ASSERT_DELTA(concs[1], 1.0,   1e-12);
            TS_ASSERT_DELTA(concs[2], 10.0,  1e-12);

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(),3u);
            TS_ASSERT_DELTA(apd90s[0], 301.4617, 1e-3);
            TS_ASSERT_DELTA(apd90s[1], 302.6703, 1e-3);
            TS_ASSERT_DELTA(apd90s[2], 311.27, 1e-2);

            TS_ASSERT_THROWS_THIS(methods.GetApd90CredibleRegions(),
                "There was no Lookup Table available for credible interval calculations with these settings.");

        }
    }

    void TestWithConfidenceIntervalOutputs(void) throw(Exception)
    {
        // Test a couple of Exceptions
        {
            CommandLineArgumentsMocker wrapper("--model 1 --plasma-concs 1 10 --pic50-ik1 4.5 --credible-intervals --plasma-conc-logscale false");
            ApPredictMethods methods;
            TS_ASSERT_THROWS_THIS(methods.Run(),
                "Lookup table (for --credible-intervals) is currently only including IKr, IKs, INa and ICaL block, you have specified additional ones so quitting.");
        }

        {
            CommandLineArgumentsMocker wrapper("--model 1 --plasma-concs 1 10 --pic50-herg 4.5 --credible-intervals --plasma-conc-logscale false");
            ApPredictMethods methods;
            TS_ASSERT_THROWS_THIS(methods.Run(),
                "No argument --pic50-spread-herg has been provided. Cannot calculate credible intervals without this.");
        }

        // Test a simple hERG block with Shannon
        {
            CommandLineArgumentsMocker wrapper("--model 1 --plasma-concs 1 10 --pic50-herg 4.5 --pic50-spread-herg 0.15 --credible-intervals --plasma-conc-logscale false");

            ApPredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();

            TS_ASSERT_EQUALS(concs.size(),3u);
            TS_ASSERT_DELTA(concs[0], 0.0,   1e-12);
            TS_ASSERT_DELTA(concs[1], 1.0,   1e-12);
            TS_ASSERT_DELTA(concs[2], 10.0,  1e-12);

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(),3u);
            TS_ASSERT_DELTA(apd90s[0], 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90s[1], 213.9762, 1e-3);
            TS_ASSERT_DELTA(apd90s[2], 229.7436, 1e-3);

            std::vector<std::pair<double,double> > apd90_credible_regions = methods.GetApd90CredibleRegions();
            TS_ASSERT_EQUALS(apd90_credible_regions.size(),3u);
            TS_ASSERT_DELTA(apd90_credible_regions[0].first, 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[0].second, 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[1].first, 212.5283, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[1].second, 219.0559, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[2].first, 217.5297, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[2].second, 259.3450, 1e-3);
        }

        // We should get much reduced credible regions with repeated pIC50 values
        {
            CommandLineArgumentsMocker wrapper("--model 1 --plasma-concs 1 10 --pic50-herg 4.5 4.4 4.6 4.5 --pic50-spread-herg 0.15 --credible-intervals --plasma-conc-logscale false");

            ApPredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();

            TS_ASSERT_EQUALS(concs.size(),3u);
            TS_ASSERT_DELTA(concs[0], 0.0,   1e-12);
            TS_ASSERT_DELTA(concs[1], 1.0,   1e-12);
            TS_ASSERT_DELTA(concs[2], 10.0,  1e-12);

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(),3u);
            TS_ASSERT_DELTA(apd90s[0], 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90s[1], 213.9730, 1e-3);
            TS_ASSERT_DELTA(apd90s[2], 229.7191, 1e-3);

            std::vector<std::pair<double,double> > apd90_credible_regions = methods.GetApd90CredibleRegions();
            TS_ASSERT_EQUALS(apd90_credible_regions.size(),3u);
            TS_ASSERT_DELTA(apd90_credible_regions[0].first, 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[0].second, 211.9333, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[1].first, 213.1849, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[1].second, 215.3668, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[2].first, 223.2300, 1e-3);
            TS_ASSERT_DELTA(apd90_credible_regions[2].second, 239.4459, 1e-3);
        }

        // This is a case Gef found that seems to get 'mean' predictions outside the credible region
        // It turned out this was an error we had introduced by taking the mean IC50 value (instead of median PIC50),
        // or even better the median of the inferred pIC50s from the credible interval pIC50 samples.
        // (which is what we are doing now).
        {
            CommandLineArgumentsMocker wrapper("--model 1 --pacing-freq 0.5 "
                    "--pic50-herg 4.86 4.45 0 4.1429 --hill-herg 3.22 1.57 1 1 "
                    "--pic50-spread-herg 0.15 --hill-spread-herg 0.21 "
                    "--credible-intervals --plasma-concs 30.0");

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
                    "--plasma-concs 30.0 ");

            ApPredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();
            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(concs.size(), apd90s.size());
            TS_ASSERT_DELTA(apd90s[0], 219.604, 1e-1);
            TS_ASSERT_DELTA(apd90s[1], 219.604, 1e-1);
            TS_ASSERT_DELTA(apd90s[2], 250.522, 1e-1);
        }
    }

    void TestWithAModelInAlternans(void) throw (Exception)
	{
        // We are trying to go so fast we get alternans, and then handle it nicely.
		CommandLineArgumentsMocker wrapper("--model 3 --plasma-concs 1 10 --pacing-freq 5 --plasma-conc-logscale false");

		ApPredictMethods methods;
		methods.Run();
		std::vector<double> concs = methods.GetConcentrations();

		TS_ASSERT_EQUALS(concs.size(),3u);
		TS_ASSERT_DELTA(concs[0], 0.0,   1e-12);
		TS_ASSERT_DELTA(concs[1], 1.0,   1e-12);
		TS_ASSERT_DELTA(concs[2], 10.0,  1e-12);

		std::vector<double> apd90s = methods.GetApd90s();
		TS_ASSERT_EQUALS(apd90s.size(),3u);
		TS_ASSERT_DELTA(apd90s[0], 137.3714, 1e-2);
		TS_ASSERT_DELTA(apd90s[1], 137.3714, 1e-2);
		TS_ASSERT_DELTA(apd90s[2], 137.3714, 1e-2);
    }

    void TestWithAModelTwoToOneStimAp(void) throw (Exception)
	{
        // We should go 'too fast' for O'Hara and see if we can get anything like sensible output.
        CommandLineArgumentsMocker wrapper("--model 6 --plasma-concs 1 --pacing-freq 5 --plasma-conc-logscale false");

		ApPredictMethods methods;
		methods.Run();
		std::vector<double> concs = methods.GetConcentrations();

		TS_ASSERT_EQUALS(concs.size(),2u);
		TS_ASSERT_DELTA(concs[0], 0.0,   1e-12);
		TS_ASSERT_DELTA(concs[1], 1.0,   1e-12);

		std::vector<double> apd90s = methods.GetApd90s();
		TS_ASSERT_EQUALS(apd90s.size(),2u);
		TS_ASSERT_DELTA(apd90s[0], 266.81, 1e-1); // Very sensitive to compiler changes, so high tolerance.
		TS_ASSERT_DELTA(apd90s[1], 266.81, 1e-1); // Very sensitive to compiler changes, so high tolerance.
    }

    void TestTroublesomeApCalculation(void) throw (Exception)
    {
        // This is a case that Gef was having some trouble with. It does look as if something odd is going
        // on, as it reports AP alternans, then tries to swap round the pacing and it has gone, at 10k paces.
        CommandLineArgumentsMocker wrapper("--model 1 --pacing-freq 0.5 --pic50-cal 3.0 --hill-cal 1 "
                "--pic50-herg 0 --hill-herg 1 "
                "--pic50-na 4.561 --hill-na 1 "
                "--plasma-concs 0 100.0");

        ApPredictMethods methods;
        methods.Run();
        std::vector<double> concs = methods.GetConcentrations();

        TS_ASSERT_EQUALS(concs.size(),3u);

        std::vector<double> apd90s = methods.GetApd90s();
        TS_ASSERT_EQUALS(apd90s.size(),3u);
        TS_ASSERT_DELTA(apd90s[0], 219.60, 1e-1); // Very sensitive to compiler changes, so high tolerance.
        TS_ASSERT_DELTA(apd90s[1], 219.60, 1e-1);
        TS_ASSERT_DELTA(apd90s[2], 217.469, 1e-1);
    }

};


#endif //_TESTAPPREDICTLONG_HPP_

#endif //_CHASTE_CVODE

