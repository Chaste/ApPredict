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
#ifdef CHASTE_CVODE

#ifndef _TESTTORSADEPREDICTLONG_HPP_
#define _TESTTORSADEPREDICTLONG_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/shared_ptr.hpp>
#include "CommandLineArgumentsMocker.hpp"
#include "TorsadePredictMethods.hpp"

class TestTorsadePredictLong : public CxxTest::TestSuite
{
public:
    /**
     * This test will wipe $CHASTE_TEST_OUTPUT/ApPredict_output/
     *
     * The tests overwrite CommandLineArguments and does some standard
     * simulations to check things are working OK...
     */
    void TestSomeFavouriteCompounds(void)
    {
        // Test a simple hERG block
        {
            // Torsade predict always uses the Grandi model, so no --model flag here.
            CommandLineArgumentsMocker wrapper("--plasma-concs 1 10 --pic50-herg 5.1 --plasma-conc-logscale false");

            TorsadePredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();

            TS_ASSERT_EQUALS(concs.size(), 3u);
            TS_ASSERT_DELTA(concs[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(concs[1], 1.0, 1e-12);
            TS_ASSERT_DELTA(concs[2], 10.0, 1e-12);

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(), 3u);
            TS_ASSERT_DELTA(apd90s[0], 286.4674, 1e-3);
            TS_ASSERT_DELTA(apd90s[1], 291.2857, 1e-3);
            TS_ASSERT_DELTA(apd90s[2], 313.4771, 1e-3);

            std::vector<unsigned> tdp_predictions = methods.GetTorsadePredictions();
            TS_ASSERT_EQUALS(tdp_predictions.size(), 3u);
            TS_ASSERT_EQUALS(tdp_predictions[0], 4u); // Cat 4
            TS_ASSERT_EQUALS(tdp_predictions[1], 3u); // Cat 3
            TS_ASSERT_EQUALS(tdp_predictions[2], 2u); // Cat 1/2
        }

        // Add a check for sensible behaviour when APD calculation fails.
        {
            CommandLineArgumentsMocker wrapper("--pic50-herg 7 --pic50-na 5 --pic50-cal 6 --hill-na 2 --hill-cal 1.5 --plasma-concs 0 100");

            TorsadePredictMethods methods;
            methods.Run();
            std::vector<double> concs = methods.GetConcentrations();

            TS_ASSERT_EQUALS(concs.size(), 3u);
            TS_ASSERT_DELTA(concs[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(concs[1], 0.001, 1e-12);
            TS_ASSERT_DELTA(concs[2], 100, 1e-12);

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(), 3u);
            TS_ASSERT_DELTA(apd90s[0], 286.4674, 1e-3);
            TS_ASSERT_DELTA(apd90s[1], 286.879, 1e-3);
            TS_ASSERT_EQUALS(std::isnan(apd90s[2]), true);

            std::vector<unsigned> tdp_predictions = methods.GetTorsadePredictions();
            TS_ASSERT_EQUALS(tdp_predictions.size(), 3u);
            TS_ASSERT_EQUALS(tdp_predictions[0], 4u); // Cat 4
            TS_ASSERT_EQUALS(tdp_predictions[1], 4u); // Cat 4
            TS_ASSERT_EQUALS(tdp_predictions[2], 2u); // Cat 1/2 - APD calc failure is high risk!
        }
    }
};

#endif //_TESTTORSADEPREDICTLONG_HPP_

#endif //_CHASTE_CVODE
