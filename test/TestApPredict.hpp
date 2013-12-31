/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef _TESTAPPREDICT_HPP_
#define _TESTAPPREDICT_HPP_

#include <cxxtest/TestSuite.h>
#include <ctime>

#include <boost/shared_ptr.hpp>
#include "CommandLineArgumentsMocker.hpp"
#include "ApPredictMethods.hpp"

class TestApPredict : public CxxTest::TestSuite
{
public:

    /**
     *
     * This test will wipe $CHASTE_TEST_OUTPUT/ApPredict_output/
     *
     * The first test overwrites CommandLineArguments and checks exceptions are thrown correctly.
     */
    void TestSomeExceptions(void) throw (Exception)
    {
        ApPredictMethods methods;
        // Check some exceptions are thrown correctly...
        {
            CommandLineArgumentsMocker wrapper("--plasma-concs 1 10 --pic50-herg 3");


            TS_ASSERT_THROWS_THIS(methods.Run(),
                    "Argument \"--model <index>\" is required");
        }
        {
            CommandLineArgumentsMocker wrapper("--model 2");

            TS_ASSERT_THROWS_THIS(methods.Run(),
                    "Argument \"--plasma-conc-high <concentration in uM>\" or \"--plasma-concs <concentrations in uM>\" is required");
        }
    }

    /**
     * This last test should emulate the standalone executable and read your command line arguments.
     */
    void TestDrugAffectByVaryingConductances(void) throw (Exception)
    {
        //////////// DEFINE PARAMETERS ///////////////
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments, and
        //char **argv = *(p_args->p_argv); // is a char** of them.
        unsigned num_args = argc-1;
        std::cout << "# " << num_args << " arguments supplied.\n" << std::flush;

        if (num_args ==0 || CommandLineArguments::Instance()->OptionExists("--help"))
        {
            std::cerr << ApPredictMethods::PrintArguments() << std::flush;
            return;
        }
        ApPredictMethods methods; // No Torsade predictions.
        methods.Run();
    }

};


#endif //_TESTAPPREDICT_HPP_

#endif //_CHASTE_CVODE

