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

#ifndef TESTMAKEALOOKUPTABLE_HPP_
#define TESTMAKEALOOKUPTABLE_HPP_

#include <boost/shared_ptr.hpp>
#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "FileFinder.hpp"
#include "LookupTableGenerator.hpp"
#include "SetupModel.hpp"
#include "SingleActionPotentialPrediction.hpp"

class TestMakeALookupTable : public CxxTest::TestSuite
{
private:
    std::string mFileName;
    std::string mArchiveFilename;

public:
    void TestGenerateAndArchiveALookupTable() throw(Exception)
    {
        // This elaborate construction contains possible tuples of:
        // * 0. command line argument and filename
        // * 1. oxmeta variable name
        // * 2. whether it is to be included in this table.
        std::vector<std::tuple<std::string, std::string, bool> > channels_included{
            std::tuple<std::string, std::string, bool>{ "hERG", "membrane_rapid_delayed_rectifier_potassium_current_conductance", false },
            std::tuple<std::string, std::string, bool>{ "IKs", "membrane_slow_delayed_rectifier_potassium_current_conductance", false },
            std::tuple<std::string, std::string, bool>{ "INa", "membrane_fast_sodium_current_conductance", false },
            std::tuple<std::string, std::string, bool>{ "ICaL", "membrane_L_type_calcium_current_conductance", false },
            std::tuple<std::string, std::string, bool>{ "Ito", "membrane_transient_outward_current_conductance", false },
            std::tuple<std::string, std::string, bool>{ "INaL", "membrane_persistent_sodium_current_conductance", false },
            std::tuple<std::string, std::string, bool>{ "IK1", "membrane_inward_rectifier_potassium_current_conductance", false }
        };

        unsigned model_index;

        if (*(CommandLineArguments::Instance()->p_argc) == 1 || CommandLineArguments::Instance()->OptionExists("--help"))
        {
            std::cout << "This will make a Please input three arguments:\n"
                         " * --model <index> (with an index as per SetupModel between 1 and 8)\n"
                         " * --hertz <freq>  (the pacing frequency in Hertz - defaults to 1Hz)\n"
                         " then a list of ion channels that you would like to block:\n"
                         " * --channels <space separated list> (choice of: hERG, ICaL, INa, IKs, Ito, INaL, IK1)\n"
                      << std::flush;
            return;
        }

        if (CommandLineArguments::Instance()->OptionExists("--model"))
        {
            model_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--model");
        }
        else
        {
            EXCEPTION("Please specify model index with --model.");
        }

        double hertz = 1.0;
        if (CommandLineArguments::Instance()->OptionExists("--hertz"))
        {
            hertz = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--hertz");
        }

        if (!CommandLineArguments::Instance()->OptionExists("--channels"))
        {
            EXCEPTION("Please enter at least 1 channel name from the list: hERG, ICaL, INa, IKs, Ito, INaL, IK1.");
        }

        // Create the model.
        SetupModel setup_model(hertz, model_index);
        const std::string model_name = setup_model.GetModel()->GetSystemName();

        // The screened ion channel should really be fast Ito in models that have it separately, rather than generic Ito.
        if (setup_model.GetModel()->HasParameter("membrane_fast_transient_outward_current_conductance"))
        {
            std::get<1>(channels_included[4]) = "membrane_fast_transient_outward_current_conductance";
        }

        std::vector<std::string> channels_requested = CommandLineArguments::Instance()->GetStringsCorrespondingToOption("--channels");
        unsigned num_channels = 0u;
        for (unsigned i = 0; i < channels_requested.size(); i++)
        {
            bool found = false;
            for (auto it = channels_included.begin(); it != channels_included.end(); ++it)
            {
                if (channels_requested[i] == std::get<0>(*it))
                {
                    std::get<2>(*it) = true;
                    num_channels++;
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                EXCEPTION("--channel option '" << channels_requested[i] << "' unrecognised. Please check spelling (see --help).");
            }
        }

        // Check the model has the desired conductances present -- quit if not.
        for (auto it = channels_included.begin(); it != channels_included.end(); ++it)
        {
            // If we requested this channel
            if (std::get<2>(*it) == true)
            {
                // If the model doesn't have this parameter, throw an error.
                if (!setup_model.GetModel()->HasParameter(std::get<1>(*it)))
                {
                    EXCEPTION("The model " << model_name << " does not have the parameter '" << std::get<1>(*it)
                                           << "' (or it is not tagged in the CellML). Cannot generate lookup table - quitting.");
                }
            }
        }

        std::stringstream file_name;
        file_name << model_name << "_";
        OutputFileHandler handler("LookupTables", false); // Don't wipe the output folder!

        file_name << num_channels << "d";

        for (auto it = channels_included.begin(); it != channels_included.end(); ++it)
        {
            if (std::get<2>(*it) == true)
            {
                file_name << "_" << std::get<0>(*it);
            }
        }
        file_name << "_" << hertz << "Hz";
        mFileName = file_name.str();

        mArchiveFilename = handler.GetOutputDirectoryFullPath() + model_name + "/" + mFileName + "_generator.arch";

        FileFinder archive_file(mArchiveFilename, RelativeTo::Absolute);

        AbstractUntemplatedLookupTableGenerator* p_generator;
        if (archive_file.IsFile())
        {
            // Create a pointer to the input archive
            std::ifstream ifs(mArchiveFilename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> p_generator;
        }
        else
        {
            // Set up a new generator
            if (num_channels == 1u)
            {
                p_generator = new LookupTableGenerator<1u>(model_index, mFileName, "LookupTables/" + model_name);
            }
            else if (num_channels == 2u)
            {
                p_generator = new LookupTableGenerator<2u>(model_index, mFileName, "LookupTables/" + model_name);
            }
            else if (num_channels == 3u)
            {
                p_generator = new LookupTableGenerator<3u>(model_index, mFileName, "LookupTables/" + model_name);
            }
            else if (num_channels == 4u)
            {
                p_generator = new LookupTableGenerator<4u>(model_index, mFileName, "LookupTables/" + model_name);
            }
            else if (num_channels == 5u)
            {
                p_generator = new LookupTableGenerator<5u>(model_index, mFileName, "LookupTables/" + model_name);
            }
            else
            {
                EXCEPTION("Something has gone wrong! Check test code.");
            }

            p_generator->SetPacingFrequency(hertz);

            for (auto it = channels_included.begin(); it != channels_included.end(); ++it)
            {
                if (std::get<2>(*it) == true)
                {
                    p_generator->SetParameterToScale(std::get<1>(*it), 0.0, 1.0);
                }
            }

            double tolerance = 1.0; //ms
            if (num_channels <= 2u)
            {
                // We can probably afford to be a bit more precise...
                tolerance = 0.5;
            }
            p_generator->AddQuantityOfInterest(Apd90, tolerance /*ms*/); // QoI and tolerance
            //p_generator->AddQuantityOfInterest(Apd50, 1 /*ms*/); // QoI and tolerance
            p_generator->SetMaxNumPaces(30u * 60u); // This is 30 minutes of 1Hz pacing
            p_generator->SetMaxVariationInRefinement(5u); // This prevents over-refining in one area.
        }

        const unsigned start_evaluations = p_generator->GetNumEvaluations();
        std::cout << "Started with " << start_evaluations << " evaluations.\n";
        const unsigned max_num_evaluations = 2000000;
        const unsigned evaluations_per_checkpoint = 500;

        for (unsigned i = start_evaluations;
             i < max_num_evaluations;
             i += evaluations_per_checkpoint)
        {
            // Run some evaluations
            p_generator->SetMaxNumEvaluations(i + evaluations_per_checkpoint);
            bool converged = p_generator->GenerateLookupTable();

            // Overwrite archive entry
            {
                AbstractUntemplatedLookupTableGenerator* const p_arch_generator = p_generator;

                std::ofstream ofs(mArchiveFilename.c_str());
                boost::archive::text_oarchive output_arch(ofs);

                output_arch << p_arch_generator;
            }

            if (converged)
            {
                break;
            }
        }
        delete p_generator;
    }
};

#endif // TESTMAKEALOOKUPTABLE_HPP_
