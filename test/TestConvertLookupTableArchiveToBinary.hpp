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

#ifndef TESTCONVERTLOOKUPTABLEARCHIVETOBINARYGENERATOR_HPP_
#define TESTCONVERTLOOKUPTABLEARCHIVETOBINARYGENERATOR_HPP_

#include <boost/shared_ptr.hpp>
#include <cxxtest/TestSuite.h>

// This is defined in the hostconfig file,
// enable it if boost_iostreams has been added to boost libraries.
#ifdef CHASTE_BOOST_IOSTREAMS
// For compressing output
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#endif

#include "CheckpointArchiveTypes.hpp"

#include "FileFinder.hpp"
#include "LookupTableGenerator.hpp"
#include "SetupModel.hpp"
#include "SingleActionPotentialPrediction.hpp"

#include "Timer.hpp"

/**
 * Here we want to load an ascii archive file, then convert it into binary for this system.
 *
 * The ascii archive that is  required may be acquired by running the ApPredict executable with
 * ./ApPredict --model 1 --credible-intervals [any drug options]
 *
 * or you can download it from
 * https://cardiac.nottingham.ac.uk/lookup_tables/
 * (will require unpacking)
 *
 * Finally we want to compare times for loading the various archive files from disk.
 *
 * If you upgrade boost, then the binary archive might become corrupt for your new version.
 * In this case, simply delete it, and this test will make it again with the new system boost.
 */
class TestConvertLookupTableArchiveToBinary : public CxxTest::TestSuite
{
private:
    boost::shared_ptr<FileFinder> mpAsciiArchiveFile;
    boost::shared_ptr<FileFinder> mpBinaryArchiveFile;
    boost::shared_ptr<FileFinder> mpCompressedAsciiArchiveFile;
    boost::shared_ptr<FileFinder> mpCompressedBinaryArchiveFile;

public:
    void TestDownloadAsciiArchiveIfMissing()
    {
        SetupModel setup(1.0, 1u); // Required to get archiving working with libraries according to #2969. Very undesirable!!

        // Just hard-code to this for now.
        std::string lookup_table_name = "shannon_wang_puglisi_weber_bers_2004_model_updated_4d_hERG_IKs_INa_ICaL_1Hz_generator";

        mpAsciiArchiveFile.reset(new FileFinder(lookup_table_name + ".arch", RelativeTo::CWD));
        mpBinaryArchiveFile.reset(new FileFinder(lookup_table_name + "_BINARY.arch", RelativeTo::CWD));
        mpCompressedAsciiArchiveFile.reset(new FileFinder(lookup_table_name + "_COMPRESSED.arch", RelativeTo::CWD));
        mpCompressedBinaryArchiveFile.reset(new FileFinder(lookup_table_name + "_BINARY_COMPRESSED.arch", RelativeTo::CWD));

        if (!mpAsciiArchiveFile->IsFile())
        {
            std::string lookup_table_URL = "https://cardiac.nottingham.ac.uk/lookup_tables/" + lookup_table_name + ".arch.tgz";
            std::cout << "\n\nAttempting to download an action potential lookup table from:\n"
                      << lookup_table_URL << "\n\n";
            EXPECT0(system, "wget " + lookup_table_URL);
            std::cout << "Download succeeded, unpacking...\n";
            EXPECT0(system, "tar xzf " + lookup_table_name + ".arch.tgz");
            std::cout << "Unpacking succeeded, removing .tgz file...\n";
            EXPECT0(system, "rm -f " + lookup_table_name + ".arch.tgz");
        }
    }

    void TestCreateTheVariousArchiveTypes()
    {
        if (!mpAsciiArchiveFile->IsFile())
        {
            WARNING("Ascii archive is not present, can't create other archives!");
            return;
        }

        // Create Binary archive if we don't already have one.
        if (!mpBinaryArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from ascii file into memory, this can take a few seconds..." << std::flush;

            // Create a pointer to the input archive
            std::ifstream ifs((mpAsciiArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            AbstractUntemplatedLookupTableGenerator* p_generator;
            input_arch >> p_generator;

            std::cout << " loaded.\nConverting to binary.\n";

            // Now archive in Binary format
            AbstractUntemplatedLookupTableGenerator* const p_arch_generator = p_generator;
            std::ofstream binary_ofs(mpBinaryArchiveFile->GetAbsolutePath().c_str(), std::ios::binary);
            boost::archive::binary_oarchive output_arch(binary_ofs);

            output_arch << p_arch_generator;
        }

#ifdef CHASTE_BOOST_IOSTREAMS
        // Now archive in compressed format?
        if (!mpCompressedAsciiArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from binary file into memory, this can take a few seconds..." << std::flush;
            Timer::Reset();

            // Create a pointer to the input archive
            std::ifstream ifs((mpBinaryArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);
            boost::archive::binary_iarchive input_arch(ifs);

            // restore from the archive
            AbstractUntemplatedLookupTableGenerator* p_generator;
            input_arch >> p_generator;

            std::cout << " loaded. It took " << Timer::GetElapsedTime() << "s.\nConverting to compressed format...\n";
            Timer::Reset();

            // Now archive in Compressed format
            AbstractUntemplatedLookupTableGenerator* const p_arch_generator = p_generator;

            // Create output file as normal.
            std::ofstream ofs(mpCompressedAsciiArchiveFile->GetAbsolutePath().c_str());

            // But instead of passing boost archive directly to this, pass boost archive
            // to a new boost filtered_ostream, and add a compressor to the filtered_ostream.
            boost::iostreams::filtering_ostream out;
            boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
            out.push(boost::iostreams::zlib_compressor(zp));
            out.push(ofs);

            // As normal, but pass the filtering_ostream instead of the ostream.
            boost::archive::text_oarchive oa(out);
            oa << p_arch_generator;

            std::cout << " done.\nIt took " << Timer::GetElapsedTime() << "s\n";
        }

        // Now archive in compressed binary format?
        if (!mpCompressedBinaryArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from binary file into memory, this can take a few seconds..." << std::flush;
            Timer::Reset();

            // Create a pointer to the input archive
            std::ifstream ifs((mpBinaryArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);
            boost::archive::binary_iarchive input_arch(ifs);

            // restore from the archive
            AbstractUntemplatedLookupTableGenerator* p_generator;
            input_arch >> p_generator;

            std::cout << " loaded. It took " << Timer::GetElapsedTime() << "s.\nConverting to compressed format...\n";
            Timer::Reset();

            // Now archive in Compressed format
            AbstractUntemplatedLookupTableGenerator* const p_arch_generator = p_generator;

            // Create output file as usual
            std::ofstream ofs(mpCompressedBinaryArchiveFile->GetAbsolutePath().c_str());

            // But instead of passing boost archive directly to this, pass boost archive
            // to a new boost filtered_ostream, and add a compressor to the filtered_ostream.
            boost::iostreams::filtering_ostream out;
            boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
            out.push(boost::iostreams::zlib_compressor(zp));
            out.push(ofs);

            // As normal, but pass the filtering_ostream instead of the ostream.
            boost::archive::binary_oarchive oa(out);
            oa << p_arch_generator;

            std::cout << " done.\nIt took " << Timer::GetElapsedTime() << "s\n";
        }
#endif // CHASTE_BOOST_IOSTREAMS
    }

    void TestTimeLoadingFromEachArchive()
    {
        if (mpAsciiArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from ascii file into memory, this can take a few seconds..." << std::flush;
            Timer::Reset();

            // Create a pointer to the input archive
            std::ifstream ifs((mpAsciiArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            AbstractUntemplatedLookupTableGenerator* p_generator;
            input_arch >> p_generator;

            std::cout << " loaded.\nIt took " << Timer::GetElapsedTime() << "s\n";
        }
        else
        {
            WARNING("Ascii archive is not present, not testing load speed!");
        }

        if (mpBinaryArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from binary file into memory, this can take a few seconds..." << std::flush;
            Timer::Reset();

            // Create a pointer to the input archive
            std::ifstream ifs((mpBinaryArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);
            boost::archive::binary_iarchive input_arch(ifs);

            // restore from the archive
            AbstractUntemplatedLookupTableGenerator* p_generator;
            input_arch >> p_generator;

            std::cout << " loaded.\nIt took " << Timer::GetElapsedTime() << "s\n";
        }
        else
        {
            WARNING("Binary archive is not present, not testing load speed!");
        }

#ifdef CHASTE_BOOST_IOSTREAMS
        if (mpCompressedAsciiArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from compressed file into memory, this can take a few seconds..." << std::flush;
            Timer::Reset();

            // Create a pointer to the input archive
            std::ifstream ifs((mpCompressedAsciiArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);

            // Set up the compressed reader
            boost::iostreams::filtering_istream in;
            boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
            in.push(boost::iostreams::zlib_decompressor(zp));
            in.push(ifs);

            // Copy compressed read to the standard boost archive input stream
            boost::archive::text_iarchive input_arch(in);

            // restore from the archive
            AbstractUntemplatedLookupTableGenerator* p_generator;
            input_arch >> p_generator;

            std::cout << " loaded.\nIt took " << Timer::GetElapsedTime() << "s\n";
        }
        else
        {
            WARNING("Compressed archive is not present, not testing load speed!");
        }

        if (mpCompressedBinaryArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from compressed binary file into memory, this can take a few seconds..." << std::flush;
            Timer::Reset();

            // Create a pointer to the input archive
            std::ifstream ifs((mpCompressedBinaryArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);

            // Set up the compressed reader
            boost::iostreams::filtering_istream in;
            boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
            in.push(boost::iostreams::zlib_decompressor(zp));
            in.push(ifs);

            // Copy compressed read to the standard boost archive input stream
            boost::archive::binary_iarchive input_arch(in);

            // restore from the archive
            AbstractUntemplatedLookupTableGenerator* p_generator;
            input_arch >> p_generator;

            std::cout << " loaded.\nIt took " << Timer::GetElapsedTime() << "s\n";
        }
        else
        {
            WARNING("Compressed binary archive is not present, not testing load speed!");
        }
#endif // CHASTE_BOOST_IOSTREAMS
    }
};

#endif // TESTCONVERTLOOKUPTABLEARCHIVETOBINARYGENERATOR_HPP_
