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

#ifndef TESTCONVERTLOOKUPTABLEARCHIVETOBINARYGENERATOR_HPP_
#define TESTCONVERTLOOKUPTABLEARCHIVETOBINARYGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

// For compressing output
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include "FileFinder.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "LookupTableGenerator.hpp"

const unsigned TABLE_DIM = 4u;

/**
 * Here we want to load an ascii archive file, then convert it into binary for this system.
 *
 * The ascii archive that is  required may be acquired by running the ApPredict executable with
 * ./ApPredict --model 1 --credible-intervals [any drug options]
 *
 * or you can download it from
 * http://www.cs.ox.ac.uk/people/gary.mirams/files/shannon_wang_puglisi_weber_bers_2004_model_updated_4d_hERG_IKs_INa_ICaL_1Hz_generator.arch.tgz
 * (will require unpacking)
 *
 * Finally we want to compare times for loading the various archive files from disk.
 */
class TestConvertLookupTableArchiveToBinary : public CxxTest::TestSuite
{
private:
    boost::shared_ptr<FileFinder> mpAsciiArchiveFile;
    boost::shared_ptr<FileFinder> mpBinaryArchiveFile;
    boost::shared_ptr<FileFinder> mpCompressedArchiveFile;

public:
    void TestConvertFileToBinary() throw (Exception)
    {
        // Just hard-code to this for now.
        mpAsciiArchiveFile.reset(new FileFinder("shannon_wang_puglisi_weber_bers_2004_model_updated_4d_hERG_IKs_INa_ICaL_1Hz_generator.arch", RelativeTo::ChasteSourceRoot));
        mpBinaryArchiveFile.reset(new FileFinder("shannon_wang_puglisi_weber_bers_2004_model_updated_4d_hERG_IKs_INa_ICaL_1Hz_generator_BINARY.arch", RelativeTo::ChasteSourceRoot));
        mpCompressedArchiveFile.reset(new FileFinder("shannon_wang_puglisi_weber_bers_2004_model_updated_4d_hERG_IKs_INa_ICaL_1Hz_generator_COMPRESSED.arch", RelativeTo::ChasteSourceRoot));

        if (!mpAsciiArchiveFile->IsFile())
        {
            WARNING("Ascii archive is not present, not creating binary file!");
            return;
        }

        // Create Binary archive if we don't already have one.
        if (!mpBinaryArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from file into memory, this can take a few seconds..." << std::flush;

            // Create a pointer to the input archive
            std::ifstream ifs((mpAsciiArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            LookupTableGenerator<TABLE_DIM>* p_generator;
            input_arch >> p_generator;

            std::cout << " loaded.\nConverting to binary.\n";

            // Now archive in Binary format
            LookupTableGenerator<TABLE_DIM>* const p_arch_generator = p_generator;
            std::ofstream binary_ofs(mpBinaryArchiveFile->GetAbsolutePath().c_str(), std::ios::binary);
            boost::archive::binary_oarchive output_arch(binary_ofs);

            output_arch << p_arch_generator;
        }

//        // Now archive in compressed format?
//        if (!mpCompressedArchiveFile->IsFile())
//        {
//            std::cout << "Loading lookup table from file into memory, this can take a few seconds..." << std::flush;
//
//            // Create a pointer to the input archive
//            std::ifstream ifs((mpAsciiArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);
//            boost::archive::text_iarchive input_arch(ifs);
//
//            // restore from the archive
//            LookupTableGenerator<TABLE_DIM>* p_generator;
//            input_arch >> p_generator;
//
//            std::cout << " loaded.\nConverting to compressed format.\n";
//
//            // Now archive in Compressed format
//            LookupTableGenerator<TABLE_DIM>* const p_arch_generator = p_generator;
//
//            std::ofstream ofs(mpCompressedArchiveFile->GetAbsolutePath().c_str());
//            boost::iostreams::filtering_ostream out;
//            boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
//            out.push(boost::iostreams::zlib_compressor(zp));
//            out.push(ofs);
//            boost::archive::binary_oarchive oa(out);
//            oa << BOOST_SERIALIZATION_NVP(p_arch_generator);
//        }

        // Now archive in compressed binary format?

    }

    void TestTimeLoadingFromEachArchive() throw (Exception)
    {
        if (mpAsciiArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from ascii file into memory, this can take a few seconds..." << std::flush;

            double start = MPI_Wtime();
            // Create a pointer to the input archive
            std::ifstream ifs((mpAsciiArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            LookupTableGenerator<TABLE_DIM>* p_generator;
            input_arch >> p_generator;

            double load_time = MPI_Wtime() - start;

            std::cout << " loaded.\nIt took " << load_time << "s\n";
        }
        else
        {
            WARNING("Ascii archive is not present, not testing load speed!");
        }

        if (mpBinaryArchiveFile->IsFile())
        {
            std::cout << "Loading lookup table from binary file into memory, this can take a few seconds..." << std::flush;

            double start = MPI_Wtime();
            // Create a pointer to the input archive
            std::ifstream ifs((mpBinaryArchiveFile->GetAbsolutePath()).c_str(), std::ios::binary);
            boost::archive::binary_iarchive input_arch(ifs);

            // restore from the archive
            LookupTableGenerator<TABLE_DIM>* p_generator;
            input_arch >> p_generator;

            double load_time = MPI_Wtime() - start;

            std::cout << " loaded.\nIt took " << load_time << "s\n";
        }
        else
        {
            WARNING("Binary archive is not present, not testing load speed!");
        }
    }
};

#endif // TESTCONVERTLOOKUPTABLEARCHIVETOBINARYGENERATOR_HPP_
