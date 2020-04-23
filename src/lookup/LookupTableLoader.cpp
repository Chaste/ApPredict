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

#include <boost/archive/archive_exception.hpp>
#include <sys/stat.h> // For system commands to download and unpack Lookup Table file.

// Chaste source includes
#include "CheckpointArchiveTypes.hpp" // Should be included first in Chaste related files.

#include "CommandLineArguments.hpp"
#include "FileFinder.hpp"

// ApPredict includes
#include "LookupTableGenerator.hpp"
#include "LookupTableLoader.hpp"

const std::string LookupTableLoader::mRemoteURL = "https://cardiac.nottingham.ac.uk/lookup_tables/";

LookupTableLoader::LookupTableLoader(const std::string& rModelName, const double& rHertz)
        : mModelName(rModelName),
          mHertz(rHertz)
{
    // Here we will attempt to use any lookup table associated with this model and
    // pacing rate.
    std::stringstream lookup_table_archive_name;
    lookup_table_archive_name << mModelName;

    // We will parse the command line arguments and work out minimal set of channels that are involved.
    // N.B. this vector should agree with the ordering in the tuple at the top of
    // TestMakeALookupTable to get the correct orderings in the filenames.
    // first -- channel name
    // second -- whether it is needed for a lookup table.
    mIdealChannelsInvolved = std::vector<std::pair<std::string, bool> >{
        std::pair<std::string, bool>{ "hERG", false },
        std::pair<std::string, bool>{ "IKs", false },
        std::pair<std::string, bool>{ "INa", false },
        std::pair<std::string, bool>{ "ICaL", false },
        std::pair<std::string, bool>{ "Ito", false },
        std::pair<std::string, bool>{ "INaL", false },
        std::pair<std::string, bool>{ "IK1", false }
    };

    // This modifies #mIdealChannelsInvolved to match command line args.
    std::string ideal_lookup_table = GetIdealTable();
    std::cout << "My ideal lookup table would be " << ideal_lookup_table << std::endl;

    LoadTableFromLocalBoostArchive(ideal_lookup_table);

    // Only continue with the logic if the local ideal table wasn't loaded.
    if (mpLookupTable == nullptr)
    {
        std::string best_lookup_table = "";

        // Get lists of available files
        std::vector<std::string> website_list = GetManifestOfTablesOnGarysWebsite();
        std::vector<std::string> local_list = GetManifestOfLocalTablesInCwd();

        // Get all plausible combos that might actually work.
        std::vector<std::string> possible_list = GenerateAllCompatibleTables();
        for (unsigned i = 0; i < possible_list.size(); i++)
        {
            //std::cout << "Possible: " << possible_list[i] << std::endl;

            // See if it appears in any local files?
            if (std::find(local_list.begin(), local_list.end(), possible_list[i]) != local_list.end())
            {
                std::cout << "Local lookup table found for " << possible_list[i] << std::endl;
                best_lookup_table = possible_list[i];
                break;
            }

            // Or web files?
            if (std::find(website_list.begin(), website_list.end(), possible_list[i]) != website_list.end())
            {
                std::cout << "Web lookup table found for " << possible_list[i] << std::endl;
                best_lookup_table = possible_list[i];
                // Download and unpack it
                DownloadAndUnpack(best_lookup_table);
                break;
            }
        }

        if (best_lookup_table != "")
        {
            LoadTableFromLocalBoostArchive(best_lookup_table);
        }
        else
        {
            WARNING("No lookup table is available, please run without --credible-intervals.");
        }
    }
}

std::string LookupTableLoader::GetIdealTable()
{
    // Parse the inputs
    CommandLineArguments* p_args = CommandLineArguments::Instance();
    unsigned ideal_dimension = 0u;

    // can't read from ApPredictMethods' mShortNames since they are in a different order, for historical reasons!!
    // so redefined here (not ideal).
    std::vector<std::string> command_line_names{ "herg", "iks", "na", "cal", "ito", "nal", "ik1" };

    // Got through command line args in order that we want them in the file name
    std::string ideal_lookup_table;
    for (unsigned i = 0; i < command_line_names.size(); i++)
    {
        if (p_args->OptionExists("--ic50-" + command_line_names[i])
            || p_args->OptionExists("--pic50-" + command_line_names[i]))
        {
            mIdealChannelsInvolved[i].second = true;
            ideal_lookup_table += "_" + mIdealChannelsInvolved[i].first;
            ideal_dimension++;
        }
    }

    // Add suffix
    std::stringstream hertz;
    hertz << mHertz;
    ideal_lookup_table = mModelName + "_" + std::to_string(ideal_dimension) + "d" + ideal_lookup_table + "_" + hertz.str() + "Hz_generator";

    return ideal_lookup_table;
}

void LookupTableLoader::LoadTableFromLocalBoostArchive(const std::string& rLookupTableBaseName)
{
    // See if there is a table available in current working directory.
    FileFinder ascii_archive_file(rLookupTableBaseName + ".arch", RelativeTo::AbsoluteOrCwd);
    FileFinder binary_archive_file(rLookupTableBaseName + "_BINARY.arch", RelativeTo::AbsoluteOrCwd);

    if (!ascii_archive_file.IsFile() && !binary_archive_file.IsFile())
    {
        // Neither local ascii nor binary are present.
        return;
    }

    // First we try loading the binary version of the archive, if it exists.
    if (binary_archive_file.IsFile())
    {
        std::cout << "Loading lookup table from binary archive into memory, this "
                     "can take a few seconds..."
                  << std::flush;
        Timer::Reset();

        // Create a pointer to the input archive
        std::ifstream ifs((binary_archive_file.GetAbsolutePath()).c_str(),
                          std::ios::binary);
        boost::archive::binary_iarchive input_arch(ifs);

        // Restore from the archive
        AbstractUntemplatedLookupTableGenerator* p_generator;
        input_arch >> p_generator;
        mpLookupTable.reset(p_generator);

        std::cout << " loaded in " << Timer::GetElapsedTime()
                  << " secs.\nLookup table is available for generation of credible "
                     "intervals.\n";

        // Since loading the binary archive works, we can try and get rid of the
        // ascii one to clean up.
        if (ascii_archive_file.IsFile())
        {
            try
            {
                // The ascii file is not in a testoutput folder so we need to over-ride
                // our usual safety checks.
                ascii_archive_file.DangerousRemove();
                std::cout << "Ascii lookup table archive file removed to tidy up, will "
                             "use the binary one in future."
                          << std::endl;
            }
            catch (Exception& e)
            {
                WARNING("Could not remove ascii lookup table archive, error was: "
                        << e.GetMessage() << "\nSimulations continued anyway.");
            }
        }

        // We have finished and loaded the generator.
        return;
    }

    // If there is no binary archive, then try loading the ascii version and
    // creating a binary one for next time.
    if (ascii_archive_file.IsFile())
    {
        std::cout << "Loading lookup table from file into memory, this can take a "
                     "few seconds..."
                  << std::flush;
        Timer::Reset();

        // Create a pointer to the input archive
        std::ifstream ifs((ascii_archive_file.GetAbsolutePath()).c_str(),
                          std::ios::binary);

        // restore from the archive
        AbstractUntemplatedLookupTableGenerator* p_generator;
        try
        {
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_generator;
        }
        catch (boost::archive::archive_exception& e)
        {
            if (std::string(e.what()) == "unsupported version")
            {
                EXCEPTION(
                    "The lookup table archive was created on a newer version of boost, "
                    "please upgrade your boost to the latest supported by this version "
                    "of Chaste.");
            }
            else
            {
                EXCEPTION("Error in loading Lookup Table from boost archive: '"
                          << e.what() << "'.");
            }
        }
        mpLookupTable.reset(p_generator);

        std::cout << " loaded in " << Timer::GetElapsedTime() << " secs."
                                                                 "\nLookup table is available for generation of credible intervals.\n";

        try
        {
            std::cout << "Saving a binary version of the archive for faster loading next time..." << std::flush;

            // Save a binary version to speed things up next time round.
            AbstractUntemplatedLookupTableGenerator* const p_arch_generator = mpLookupTable.get();

            std::ofstream binary_ofs(binary_archive_file.GetAbsolutePath().c_str(),
                                     std::ios::binary);
            boost::archive::binary_oarchive output_arch(binary_ofs);

            output_arch << p_arch_generator;

            std::cout << "done!\n";
        }
        catch (Exception& e)
        {
            WARNING("Did not manage to create binary lookup table archive. Error was: "
                    << e.GetMessage() << "\nContinuing to use ascii archive.");
        }
    }
}

std::vector<std::string> LookupTableLoader::GetManifestOfTablesOnGarysWebsite()
{
    const std::string manifest_filename = "appredict_lookup_table_manifest.txt";
    std::vector<std::string> available_tables;

    // If no archive exists, try to download and unpack one.
    std::string mainfest_URL = mRemoteURL + manifest_filename;
    FileFinder manifest(manifest_filename, RelativeTo::AbsoluteOrCwd);

    // First check to see whether the remote manifest is accessible
    // The "--server-response" argument prints the server status and the "--spider" means don't download anything.
    std::string command = "wget --server-response --spider " + mainfest_URL;
    int return_code = system(command.c_str());
    if (return_code != 0)
    {
        std::cout << "Could not find the remote manifest of available Lookup Tables on the web, "
                     "we either don't have web access or the lookup table host server is down..."
                  << std::endl;
        return available_tables;
    }

    try
    {
        if (manifest.IsFile())
        {
            std::cout << "\n\nAttempting to overwrite local lookup table manifest with the latest from:\n"
                      << mainfest_URL << "\n\n";
        }
        else
        {
            std::cout << "\n\nAttempting to download lookup table manifest from:\n"
                      << mainfest_URL << "\n\n";
        }

        EXPECT0(system, "wget --dns-timeout=10 --connect-timeout=10 -O " + manifest_filename + " " + mainfest_URL);

        std::cout << "Download succeeded.\n";
    }
    catch (Exception& e)
    {
        std::cout << "Could not download and unpack the Lookup Table manifest, "
                     "we either don't have web access or the lookup table host server is down..."
                  << std::endl;
        return available_tables;
    }

    EXCEPT_IF_NOT(manifest.IsFile());

    std::ifstream manifest_stream(manifest_filename);
    std::string filename;

    // Check for suitable Hertz
    std::stringstream hertz;
    hertz << mHertz;
    std::string hertz_string = hertz.str();

    while (manifest_stream >> filename)
    {
        // If the manifest file starts with the model name
        if (filename.compare(0, mModelName.length(), mModelName) == 0)
        {
            // Look for 'XHz_generator'
            size_t hertz_pos = filename.find(hertz_string + "Hz_generator");
            if (hertz_pos != std::string::npos)
            {
                // Take off .arch.tgz
                size_t extension_pos = filename.find(".arch.tgz");
                available_tables.push_back(filename.substr(0, extension_pos));
            }
        }
    }
    return available_tables;
}

std::vector<std::string> LookupTableLoader::GetManifestOfLocalTablesInCwd()
{
    FileFinder cwd("", RelativeTo::AbsoluteOrCwd);

    std::vector<FileFinder> matching_files = cwd.FindMatches("*.arch");

    const std::string binary_ending = "_BINARY";

    // Check for suitable Hertz
    std::stringstream hertz;
    hertz << mHertz;
    std::string hertz_string = hertz.str();

    std::vector<std::string> available_tables;
    for (unsigned i = 0; i < matching_files.size(); i++)
    {
        std::string base_file_name = matching_files[i].GetLeafNameNoExtension();

        // If matching file ends in _BINARY
        if (0 == base_file_name.compare(base_file_name.length() - binary_ending.length(), binary_ending.length(), binary_ending))
        {
            // Strip off the "_BINARY" extension
            base_file_name = base_file_name.substr(0, base_file_name.length() - binary_ending.length());
        }

        // If the manifest file starts with the model name
        if (base_file_name.compare(0, mModelName.length(), mModelName) == 0)
        {
            // Look for 'XHz_generator'
            size_t hertz_pos = base_file_name.find(hertz_string + "Hz_generator");
            if (hertz_pos != std::string::npos)
            {
                // Only add if it isn't already present
                if (std::find(available_tables.begin(), available_tables.end(), base_file_name) == available_tables.end())
                {
                    available_tables.push_back(base_file_name);
                }
            }
        }
    }
    return available_tables;
}

std::vector<std::string> LookupTableLoader::GenerateAllCompatibleTables()
{

    // Work out how many we start with, ideal number
    unsigned ideal_num = 0;
    std::vector<bool> ideal(7, false);
    for (unsigned i = 0; i < mIdealChannelsInvolved.size(); i++)
    {
        if (mIdealChannelsInvolved[i].second == true)
        {
            ideal_num++;
            ideal[i] = true;
        }
    }

    unsigned max_num_to_add = mIdealChannelsInvolved.size() - ideal_num;

    /*
     * Now follows some fairly horrible logic to do combinations that would work.
     *
     * Probably a far simpler way to do this, but it works!
     */
    std::vector<std::vector<bool> > compatible_channels;
    compatible_channels.push_back(ideal);
    std::vector<std::vector<bool> > compatible_channels_added_this_dim;
    std::vector<std::vector<bool> > start_points_this_dimension;

    start_points_this_dimension.push_back(ideal);
    // Loop over the dimensions (add 1,2,3,4...)
    for (unsigned dim = 0; dim < max_num_to_add; dim++)
    {
        compatible_channels_added_this_dim.clear();
        for (unsigned s = 0; s < start_points_this_dimension.size(); s++)
        {
            std::vector<bool> start_point = start_points_this_dimension[s];
            // Look over the channels
            for (unsigned i = 0; i < mIdealChannelsInvolved.size(); i++)
            {
                if (mIdealChannelsInvolved[i].second == false)
                {
                    std::vector<bool> possible = start_point;
                    possible[i] = true;
                    if (std::find(compatible_channels.begin(), compatible_channels.end(), possible) == compatible_channels.end())
                    {
                        compatible_channels.push_back(possible);
                        compatible_channels_added_this_dim.push_back(possible);
                    }
                }
            }
        }
        start_points_this_dimension.clear();
        start_points_this_dimension = compatible_channels_added_this_dim;
    }

    // Convert bools to filenames
    std::vector<std::string> compatible_tables;
    for (unsigned i = 0; i < compatible_channels.size(); i++)
    {
        std::vector<bool> channels_this_combo = compatible_channels[i];
        std::string channels_text;
        unsigned dimension = 0;
        for (unsigned j = 0; j < mIdealChannelsInvolved.size(); j++)
        {
            if (channels_this_combo[j])
            {
                channels_text += "_" + mIdealChannelsInvolved[j].first;
                dimension++;
            }
        }

        std::stringstream name;
        name << mModelName << "_" << dimension << "d" << channels_text << "_" << mHertz << "Hz_generator";
        compatible_tables.push_back(name.str());
    }
    return compatible_tables;
}

void LookupTableLoader::DownloadAndUnpack(const std::string& rArchiveFileBaseName)
{
    std::string lookup_table_URL = mRemoteURL + rArchiveFileBaseName + ".arch.tgz";
    try
    {
        std::cout << "\n\nAttempting to download an action potential lookup table from:\n"
                  << lookup_table_URL << "\n\n";

        EXPECT0(system, "wget --dns-timeout=10 --connect-timeout=10 " + lookup_table_URL);

        std::cout << "Download succeeded, unpacking...\n";

        EXPECT0(system, "tar xzf " + rArchiveFileBaseName + ".arch.tgz");

        std::cout << "Unpacking succeeded, removing .tgz file...\n";

        EXPECT0(system, "rm -f " + rArchiveFileBaseName + ".arch.tgz");
    }
    catch (Exception& e)
    {
        std::cout << "Could not download and unpack the Lookup Table archive, continuing without it..." << std::endl;
    }
}

bool LookupTableLoader::IsLookupTableAvailable()
{
    return (mpLookupTable != nullptr);
}

boost::shared_ptr<AbstractUntemplatedLookupTableGenerator> LookupTableLoader::GetLookupTable()
{
    if (mpLookupTable == nullptr)
    {
        EXCEPTION("A lookup table could not be loaded.");
    }

    return mpLookupTable;
}
