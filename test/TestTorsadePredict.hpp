/*

Copyright (c) 2005-2020, University of Oxford.
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

#ifndef _TESTTORSADEPREDICT_HPP_
#define _TESTTORSADEPREDICT_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/shared_ptr.hpp>
#include "CommandLineArguments.hpp"
#include "TorsadePredictMethods.hpp"

class TestTorsadePredict : public CxxTest::TestSuite
{
public:
    void TestGetTorsadePredictions(void)
    {
        // Common sense checks - this one should be for high APD -> danger
        {
            TorsadePredictMethods ap_methods;

            TS_ASSERT_THROWS_THIS(ap_methods.MakeTorsadePredictions(),
                                  "APDs do not appear to have been recorded.");

            ap_methods.mApd90s.push_back(282.493);
            ap_methods.mApd90s.push_back(290);
            ap_methods.mApd90s.push_back(312);
            ap_methods.mApd90s.push_back(333);

            ap_methods.MakeTorsadePredictions();

            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions.size(), 4u);
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[0], 4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[1], 3u); // Cat 3
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[2], 2u); // Cat 1/2
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[3], 2u); // Cat 1/2

            // Check the get method works OK.
            TS_ASSERT_THROWS_THIS(ap_methods.GetTorsadePredictions(),
                                  "Simulation has not been run - check arguments.");

            // Hack it to make it look as if it did run properly
            ap_methods.mComplete = true;

            std::vector<unsigned> tdp_predicitons = ap_methods.GetTorsadePredictions();
            TS_ASSERT_EQUALS(tdp_predicitons.size(), ap_methods.mTorsadePredictions.size());
            for (unsigned i = 0; i < tdp_predicitons.size(); i++)
            {
                TS_ASSERT_EQUALS(tdp_predicitons[i], ap_methods.mTorsadePredictions[i]);
            }
        }

        // Common sense checks - this one should be for low APD -> safer
        {
            TorsadePredictMethods ap_methods;
            ap_methods.mApd90s.push_back(282.493);
            ap_methods.mApd90s.push_back(280);
            ap_methods.mApd90s.push_back(275);
            ap_methods.mApd90s.push_back(270);

            ap_methods.MakeTorsadePredictions();

            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions.size(), 4u);
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[0], 4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[1], 4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[2], 4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[3], 4u); // Cat 4
        }

        // Common sense checks - this one should be for very low APD -> very safe
        {
            TorsadePredictMethods ap_methods;
            ap_methods.mApd90s.push_back(282.493);
            ap_methods.mApd90s.push_back(260);
            ap_methods.mApd90s.push_back(235);
            ap_methods.mApd90s.push_back(230);

            ap_methods.MakeTorsadePredictions();

            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions.size(), 4u);
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[0], 4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[1], 5u); // Cat 5
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[2], 5u); // Cat 5
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[3], 5u); // Cat 5
        }

        // Common sense checks - this one goes very high and then low again should be dodgy all the way up.
        {
            TorsadePredictMethods ap_methods;
            ap_methods.mApd90s.push_back(282.493);
            ap_methods.mApd90s.push_back(330);
            ap_methods.mApd90s.push_back(300);
            ap_methods.mApd90s.push_back(282.493);
            ap_methods.mApd90s.push_back(230);

            ap_methods.MakeTorsadePredictions();

            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions.size(), 5u);
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[0], 4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[1], 2u); // Cat 1/2
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[2], 2u); // Cat 1/2
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[3], 2u); // Cat 1/2
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[4], 2u); // Cat 1/2
        }
    }

    /**
     *
     * This test will wipe $CHASTE_TEST_OUTPUT/TorsadePredict_output/
     *
     * Parameters can be defined at the top of this Test
     */
    void TestDrugAffectByVaryingConductances(void)
    {
        //////////// DEFINE PARAMETERS ///////////////
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments, and
        //char **argv = *(p_args->p_argv); // is a char** of them.
        unsigned num_args = argc - 1;
        std::cout << "# " << num_args << " arguments supplied.\n"
                  << std::flush;

        if (num_args == 0 || CommandLineArguments::Instance()->OptionExists("--help"))
        {
            std::cerr << TorsadePredictMethods::PrintArguments() << std::flush;
            return;
        }
        TorsadePredictMethods methods;
        methods.Run();
    }
};

#endif //_TESTTORSADEPREDICT_HPP_

#endif //_CHASTE_CVODE
