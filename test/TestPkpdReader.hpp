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

#ifndef TESTPKPDREADER_HPP_
#define TESTPKPDREADER_HPP_

#include <cxxtest/TestSuite.h>

#include "PkpdDataStructure.hpp"

/**
 * A test that checks we are reading in concentration information correctly.
 */
class TestPkpdReader : public CxxTest::TestSuite
{
public:
    void TestPkPdDataReader()
    {
        FileFinder pkpd_data_file("projects/ApPredict/test/data/pkpd_data.txt",
                                  RelativeTo::ChasteSourceRoot);
        TS_ASSERT_EQUALS(pkpd_data_file.IsFile(), true);

        PkpdDataStructure pkpd_data(pkpd_data_file);

        std::vector<std::string> times = pkpd_data.GetTimes();
        TS_ASSERT_EQUALS(times.size(), 749u);

        TS_ASSERT_EQUALS(pkpd_data.GetNumberOfPatients(), 57u);

        std::vector<double> concs = pkpd_data.GetConcentrationsForPatient(0u);
        TS_ASSERT_EQUALS(concs.size(), times.size());

        TS_ASSERT_THROWS_THIS(pkpd_data.GetConcentrationsForPatient(57u),
                              "Patient index 57 requested but there are only 57 in the data file.");

        // Test concentration interpolation methods
        TS_ASSERT_EQUALS(times[0], "0");
        TS_ASSERT_EQUALS(times.back(), "748");
        TS_ASSERT_DELTA(concs[0], 0, 1e-9);
        TS_ASSERT_DELTA(concs.back(), 0.0195722, 1e-9); // uM

        std::vector<double> concs2 = pkpd_data.GetConcentrationsForPatient(56u);
        TS_ASSERT_EQUALS(concs2.size(), times.size());

        TS_ASSERT_DELTA(concs2[0], 0, 1e-9);
        TS_ASSERT_DELTA(concs2.back(), 1.11562, 1e-9); // uM

        TS_ASSERT_DELTA(pkpd_data.GetMaximumConcentration(), 4.1515, 1e-4); // uM
    }

    void TestPkPdDataReaderDos()
    {
        FileFinder pkpd_data_file("projects/ApPredict/test/data/DosTestFile.txt",
                                  RelativeTo::ChasteSourceRoot);
        TS_ASSERT_EQUALS(pkpd_data_file.IsFile(), true);

        PkpdDataStructure pkpd_data(pkpd_data_file);

        std::vector<std::string> times = pkpd_data.GetTimes();
        TS_ASSERT_EQUALS(times.size(), 10u);

        TS_ASSERT_EQUALS(pkpd_data.GetNumberOfPatients(), 4u);
    }

    void TestPkPdDataReaderUnix()
    {
        FileFinder pkpd_data_file("projects/ApPredict/test/data/UnixTestFile.txt",
                                  RelativeTo::ChasteSourceRoot);
        TS_ASSERT_EQUALS(pkpd_data_file.IsFile(), true);

        PkpdDataStructure pkpd_data(pkpd_data_file);

        std::vector<std::string> times = pkpd_data.GetTimes();
        TS_ASSERT_EQUALS(times.size(), 10u);

        TS_ASSERT_EQUALS(pkpd_data.GetNumberOfPatients(), 4u);
    }
};

#endif // TESTPKPDREADER_HPP_
