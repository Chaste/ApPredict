/*

Copyright (c) 2005-2022, University of Oxford.
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
#include <vector>

#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"

#include "CommandLineArguments.hpp"
#include "RunHillFunctionMinimization.hpp"

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StandardStartup(&argc, &argv);

    int exit_code = ExecutableSupport::EXIT_OK;

    // You should put all the main code within a try-catch, to ensure that
    // you clean up PETSc before quitting.
    try
    {
        if (PetscTools::AmMaster())
        {
            unsigned num_parameters_to_fit = 1; // is default

            if (CommandLineArguments::Instance()->OptionExists("--num-params"))
            {
                unsigned choice = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--num-params");
                if (choice == 1u || choice == 2u )
                {
                    num_parameters_to_fit = choice;
                }
                else
                {
                    EXCEPTION("Can only fit 1 or 2 parameters, not " << choice << ".");
                    exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
                }
            }

            if (!(CommandLineArguments::Instance()->OptionExists("--concs") && CommandLineArguments::Instance()->OptionExists("--responses")))
            {
                EXCEPTION("You must supply the options:\n"
                          "*  --concs  (followed by a list of concentrations in uM).\n"
                          "*  --responses  (followed by a list of percent inhibitions.\n"
                          "optionally:\n"
                          "*  --num-params  1/2 Whether to fit one or two parameters, IC50 or IC50 & Hill.\n"
                          "(defaults to just IC50)\n"
                          "*  --hill-limits <x> <y> The minimum and maximum limits to use when fitting a Hill coefficient.\n"
                          );
                exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
            }

            std::vector<double> concentrations = CommandLineArguments::Instance()->GetDoublesCorrespondingToOption("--concs");
            std::vector<double> inhibitions = CommandLineArguments::Instance()->GetDoublesCorrespondingToOption("--responses");

            if (concentrations.size()!=inhibitions.size())
            {
                EXCEPTION("The list of concentrations and responses must be the same length,\n"
                         "they appear to be concs.size() = " << concentrations.size() << " and inhibitions.size() = " << inhibitions.size());
                exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
            }

            if (concentrations.size()==1)
            {
                // Don't try and fit two parameters to one data point.
                num_parameters_to_fit = 1u;
            }

            std::cout<< "Fit is using " << concentrations.size() << " dose-response points.\n";

            RunHillFunctionMinimization compound(concentrations,inhibitions,num_parameters_to_fit);

            if (CommandLineArguments::Instance()->OptionExists("--hill-limits"))
            {
                std::vector<double> limits = CommandLineArguments::Instance()->GetDoublesCorrespondingToOption("--hill-limits");
                if (limits.size() != 2u || limits[1] <= limits[0])
                {
                    EXCEPTION("The command line option \"--hill-limits\" must be followed by two numeric values for min and max.");
                }
                compound.SetHillLimits(limits[0], limits[1]);
            }
            std::vector<double> parameters = compound.Run();

            //Prints the calculated pIC50 value and hill coefficient to screen at end of optimisation process
            if (parameters.size()==1)
            {
                parameters.push_back(1.0);
                std::cout<<"The IC50 is " << parameters[0] << "uM, [pIC50 is "<<-log10((1e-6)*(parameters[0]))<<" (log M)]\nand the hill coefficient is 1.\n";
            }
            else
            {
                std::cout<<"The IC50 is " << parameters[0] << "uM, [pIC50 is "<<-log10((1e-6)*(parameters[0]))<<" (log M)]\nand the hill coefficient is "<<parameters[1]<<".\n";
            }
        }
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

    // Optional - write the machine info to file.
    //ExecutableSupport::WriteMachineInfoFile("machine_info");

    // End by finalizing PETSc, and returning a suitable exit code.
    // 0 means 'no error'
    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
