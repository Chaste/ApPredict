/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTACTIONPOTENTIALDOWNSAMPLER_HPP_
#define TESTACTIONPOTENTIALDOWNSAMPLER_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractDataReader.hpp"
#include "ActionPotentialDownsampler.hpp"
#include "CommandLineArgumentsMocker.hpp"
#include "FileComparison.hpp"
#include "FileFinder.hpp"

class TestActionPotentialDownsampler : public CxxTest::TestSuite
{
public:
    void TestAgainstStoredFile()
    {
        std::string output_folder = "TestActionPotentialDownsampler";
        std::string output_filename = "sample_action_potential.txt";

        std::vector<double> times;
        std::vector<double> voltages;
        // Get a full action potential trace (large file) out of one stored in the repository...
        {
            std::ifstream indata; // indata is like cin
            indata.open("projects/ApPredict/test/data/full_voltage_trace.dat"); // opens the file
            assert(indata.good());

            bool top_line_read = false; // There is a row of header info we need to skip.
            while (indata.good())
            {
                std::string this_line;
                getline(indata, this_line);

                if (!top_line_read)
                {
                    top_line_read = true;
                    continue;
                }

                if (this_line == "" || this_line == "\r")
                {
                    if (indata.eof())
                    { // If the blank line is the last line carry on OK.
                        break;
                    }
                }
                // Put the data in a stringstream and then into the vectors.
                std::stringstream line(this_line);
                double time;
                double voltage;
                line >> time;
                line >> voltage;
                times.push_back(time);
                voltages.push_back(voltage);
            }
        }

        // Options we would like to pass to the downsampler.
        double window = 500; //ms
        double stim_time = 5; //ms

        // Test with downsampling switched on (default)
        {
            ActionPotentialDownsampler sampler(output_folder, output_filename, times, voltages, window, stim_time);

            FileFinder generated_file("TestActionPotentialDownsampler/sample_action_potential.txt", RelativeTo::ChasteTestOutput);
            FileFinder reference_file("projects/ApPredict/test/data/reduced_voltage_trace.dat", RelativeTo::ChasteSourceRoot);
            TS_ASSERT(generated_file.IsFile());
            TS_ASSERT(reference_file.IsFile());

            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with downsampling switched off
        CommandLineArgumentsMocker wrapper("--no-downsampling");
        output_filename = "sample_action_potential_no_downsampling.txt";

        // These are altered just so we can directly compare with the full_voltage_trace of input data.
        stim_time = 0; //ms
        window = 2000; //ms
        {
            ActionPotentialDownsampler sampler(output_folder, output_filename, times, voltages, window, stim_time);

            FileFinder generated_file("TestActionPotentialDownsampler/sample_action_potential_no_downsampling.txt", RelativeTo::ChasteTestOutput);
            FileFinder reference_file("projects/ApPredict/test/data/full_voltage_trace.dat", RelativeTo::ChasteSourceRoot);
            TS_ASSERT(generated_file.IsFile());
            TS_ASSERT(reference_file.IsFile());

            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif // TESTACTIONPOTENTIALDOWNSAMPLER_HPP_
