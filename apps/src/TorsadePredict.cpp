/*

Copyright (c) 2005-2019, University of Oxford.
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

#include <iostream>
#include <string>

#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "CommandLineArguments.hpp"
#include "TorsadePredictMethods.hpp"

#include "PetscTools.hpp"
#include "PetscException.hpp"

/**
 * The Torsade prediction program that Mirams et al. 2011 suggests.
 *
 * Requires the file "drug_data.dat" to be in the current working directory.
 */
int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StandardStartup(&argc, &argv);
    ExecutableSupport::SetOutputDirectory("TorsadePredict_output");

    int exit_code = ExecutableSupport::EXIT_OK;

    // You should put all the main code within a try-catch, to ensure that
    // you clean up PETSc before quitting.
    try
    {
        //////////// DEFINE PARAMETERS ///////////////
        unsigned num_args = argc-1;
        std::cout << "# " << num_args << " arguments supplied.\n" << std::flush;

        if (num_args == 0 || CommandLineArguments::Instance()->OptionExists("--help"))
        {
            ExecutableSupport::PrintError(TorsadePredictMethods::PrintArguments());
            ExecutableSupport::FinalizePetsc();
            return ExecutableSupport::EXIT_BAD_ARGUMENTS;
        }
        TorsadePredictMethods methods;  // Use Grandi model, record APD90 and do Torsade Predictions.
        methods.Run();
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

    ExecutableSupport::WriteMachineInfoFile("machine_info");

    // End by finalising PETSc, and returning a suitable exit code.
    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
