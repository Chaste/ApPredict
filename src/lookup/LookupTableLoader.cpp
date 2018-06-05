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

#include <boost/archive/archive_exception.hpp>
#include <sys/stat.h> // For system commands to download and unpack Lookup Table file.

// Chaste source includes
#include "CheckpointArchiveTypes.hpp" // Should be included first in Chaste related files.

#include "CommandLineArguments.hpp"
#include "FileFinder.hpp"

// ApPredict includes
#include "LookupTableGenerator.hpp"
#include "LookupTableLoader.hpp"

LookupTableLoader::LookupTableLoader(const std::string& rModelName, const double& rHertz)
        : mWeHaveWebAccess(true),
          mModelName(rModelName),
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
        std::pair<std::string, bool>{ "INa", false },
        std::pair<std::string, bool>{ "IKs", false },
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
    if (!mpLookupTable)
    {
        std::string best_lookup_table;

        // Get list of remote files
        std::vector<std::string> website_list = GetManifestOfTablesOnGarysWebsite();

        for (unsigned i = 0; i < website_list.size(); i++)
        {
            std::cout << "On web: " << website_list[i] << std::endl;
        }
        assert(0);

        LoadTableFromLocalBoostArchive(best_lookup_table);
    }
}

std::string LookupTableLoader::GetIdealTable()
{
    // Parse the inputs
    CommandLineArguments* p_args = CommandLineArguments::Instance();
    unsigned ideal_dimension = 0u;

    // can't read from ApPredictMethods' mShortNames since they are in a different order, for historical reasons!!
    // so redefined here (not ideal).
    std::vector<std::string> command_line_names{ "herg", "na", "iks", "cal", "ito", "nal", "ik1" };

    // Got through command line args in order that we want them in the file name
    std::string ideal_lookup_table = mModelName + "_" + std::to_string(ideal_dimension) + "d";
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
    ideal_lookup_table += "_" + hertz.str() + "Hz_generator";

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

    // If no archive exists, try to download and unpack one.
    std::string mainfest_URL = "http://www.cs.ox.ac.uk/people/gary.mirams/files/" + manifest_filename;
    FileFinder manifest(manifest_filename, RelativeTo::AbsoluteOrCwd);
    try
    {
        // Clean up if this is running afresh, otherwise we use a stale archive
        // because wget copies a new one to a different place with a suffix of ".<number>"
        if (manifest.IsFile())
        {
            EXPECT0(system, "rm " + manifest_filename);
        }

        std::cout << "\n\nAttempting to download lookup table manifest from:\n"
                  << mainfest_URL << "\n\n";

        EXPECT0(system, "wget --dns-timeout=10 --connect-timeout=10 " + mainfest_URL);

        std::cout << "Download succeeded.\n";
    }
    catch (Exception& e)
    {
        std::cout << "Could not download and unpack the Lookup Table manifest, "
                     "we either don't have web access or www.cs.ox.ac.uk is down..."
                  << std::endl;
        mWeHaveWebAccess = false;
    }

    EXCEPT_IF_NOT(manifest.IsFile());

    std::ifstream manifest_stream(manifest_filename);
    std::string filename;

    // Check for suitable Hertz
    std::stringstream hertz;
    hertz << mHertz;
    std::string hertz_string = hertz.str();

    std::vector<std::string> available_tables;
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

    std::vector<FileFinder> matching_files = cwd.FindMatches(mModelName + "*_generator*.arch");

    std::vector<std::string> available_tables;
    for (unsigned i = 0; i < matching_files.size(); i++)
    {
        available_tables.push_back(matching_files[i].GetLeafNameNoExtension());
        std::cout << "Found local file " << matching_files[i].GetLeafNameNoExtension() << std::endl;
    }
    return available_tables;
}

std::vector<std::string> LookupTableLoader::GenerateAllCompatibleTables()
{
    std::vector<std::string> compatible_tables;
    return compatible_tables;
}

void LookupTableLoader::DownloadAndUnpack(const std::string& rArchiveFileBaseName)
{
    assert(mWeHaveWebAccess);

    std::string lookup_table_URL = "http://www.cs.ox.ac.uk/people/gary.mirams/files/" + rArchiveFileBaseName + ".arch.tgz";
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

boost::shared_ptr<AbstractUntemplatedLookupTableGenerator> LookupTableLoader::GetLookupTable()
{
    if (!mpLookupTable)
    {
        EXCEPTION("A lookup table could not be loaded.");
    }

    return mpLookupTable;
}
