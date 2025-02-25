/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTPKPDINTERPOLATIONS_HPP_
#define TESTPKPDINTERPOLATIONS_HPP_

#include <cxxtest/TestSuite.h>

#include "ApPredictMethods.hpp"
#include "CommandLineArgumentsMocker.hpp"
#include "NumericFileComparison.hpp"

/**
 * A test that checks we can do interpolations with ApPredict runs based on a
 * PKPD file.
 */
class TestPkpdInterpolations : public CxxTest::TestSuite
{
public:
    void TestPkpdExceptions()
    {
        CommandLineArgumentsMocker wrapper("--pkpd-file nonsense.txt --model 1");

        ApPredictMethods bad_pkpd_runner;
        TS_ASSERT_THROWS_CONTAINS(bad_pkpd_runner.Run(),
                                  "does not exist. Please give a relative or absolute path.");
    }

    void TestPkpdSimulations()
    {
        CommandLineArgumentsMocker wrapper(
            "--pkpd-file projects/ApPredict/test/data/pkpd_data.txt --model 2 --pic50-herg 6");

        ApPredictMethods pkpd_runner;

        { // Test some interpolation methods
            std::vector<double> x{ 0.0, 1.0, 2.0, 3.0 };
            std::vector<double> y{ 1.0, 1.1, -0.1, 0.0 };
            TS_ASSERT_DELTA(pkpd_runner.DoLinearInterpolation(-0.1, x, y), 1.0, 1e-6);
            TS_ASSERT_DELTA(pkpd_runner.DoLinearInterpolation(0.0, x, y), 1.0, 1e-6);
            TS_ASSERT_DELTA(pkpd_runner.DoLinearInterpolation(1.0, x, y), 1.1, 1e-6);
            TS_ASSERT_DELTA(pkpd_runner.DoLinearInterpolation(0.5, x, y), 1.05, 1e-6);
            TS_ASSERT_DELTA(pkpd_runner.DoLinearInterpolation(1.5, x, y), 0.5, 1e-6);
            TS_ASSERT_DELTA(pkpd_runner.DoLinearInterpolation(2.5, x, y), -0.05,
                            1e-6);
            TS_ASSERT_DELTA(pkpd_runner.DoLinearInterpolation(3.0, x, y), 0.0, 1e-6);
            TS_ASSERT_DELTA(pkpd_runner.DoLinearInterpolation(3.1, x, y), 0.0, 1e-6);
        }

        // Run a real simulation
        pkpd_runner.Run();

        FileFinder pkpd_results_file("ApPredict_output/pkpd_results.txt", RelativeTo::ChasteTestOutput);
        TS_ASSERT_EQUALS(pkpd_results_file.IsFile(), true);

        FileFinder pkpd_reference_file("projects/ApPredict/test/data/pkpd_results.txt", RelativeTo::ChasteSourceRoot);

        NumericFileComparison comparison(pkpd_results_file, pkpd_reference_file);
        comparison.CompareFiles(1e-2);
    }
};

#endif // TESTPKPDINTERPOLATIONS_HPP_
