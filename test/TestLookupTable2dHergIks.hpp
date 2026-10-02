/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef TESTLOOKUPTABLE2DHERGIKS_HPP_
#define TESTLOOKUPTABLE2DHERGIKS_HPP_

#include <cxxtest/TestSuite.h>
#include <iomanip>

#include "CheckpointArchiveTypes.hpp"

#include "LookupTableGenerator.hpp"
#include "OutputFileHandler.hpp"

/**
 * Generates a 2D APD90 lookup table for hERG and IKs block (ten Tusscher 2006 epi at 1Hz)
 * refined to a 1ms error tolerance, like the example plotted at the end of
 * https://github.com/Chaste/trac_archive/blob/master/issues/2366.md
 *
 * It writes the evaluated points (2d_hERG_IKs_1Hz.dat) and the table interpolated onto a
 * regular grid (2d_hERG_IKs_1Hz_interpolated.dat), so the interpolated surface and the
 * distribution of evaluated points can both be plotted.
 */
class TestLookupTable2dHergIks : public CxxTest::TestSuite
{
public:
    void TestGenerate2dHergIksTable()
    {
        const std::string output_folder = "TestLookupTable2dHergIks";
        const std::string file_name = "2d_hERG_IKs_1Hz";
        OutputFileHandler handler(output_folder); // Wipe the folder for a fresh table each time.

        LookupTableGenerator<2> generator(2u, file_name, output_folder); // Ten Tusscher 2006 epi
        generator.SetParameterToScale("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
        generator.SetParameterToScale("membrane_slow_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
        generator.AddQuantityOfInterest(Apd90, 1.0 /*ms*/); // QoI and tolerance
        generator.SetMaxNumPaces(30u * 60u); // 30 minutes of 1Hz pacing, as in TestMakeALookupTable
        generator.SetMaxVariationInRefinement(5u);
        generator.SetMaxNumEvaluations(20000u);

        bool converged = generator.GenerateLookupTable();
        TS_ASSERT(converged);
        std::cout << "Lookup table used " << generator.GetNumEvaluations() << " evaluations.\n";

        // Interpolate the table onto a regular grid for plotting the surface.
        const unsigned num_grid_points = 101u;
        std::vector<std::vector<double> > grid_points;
        for (unsigned i = 0; i < num_grid_points; i++)
        {
            for (unsigned j = 0; j < num_grid_points; j++)
            {
                grid_points.push_back(std::vector<double>{ (double)(i) / (num_grid_points - 1u),
                                                           (double)(j) / (num_grid_points - 1u) });
            }
        }
        std::vector<std::vector<double> > interpolated = generator.Interpolate(grid_points);

        out_stream p_file = handler.OpenOutputFile(file_name + "_interpolated.dat");
        *p_file << std::setprecision(8);
        *p_file << "hERG\tIKs\tAPD90\n";
        for (unsigned i = 0; i < grid_points.size(); i++)
        {
            *p_file << grid_points[i][0] << "\t" << grid_points[i][1] << "\t" << interpolated[i][0] << "\n";
        }
        p_file->close();
    }
};

#endif // TESTLOOKUPTABLE2DHERGIKS_HPP_
