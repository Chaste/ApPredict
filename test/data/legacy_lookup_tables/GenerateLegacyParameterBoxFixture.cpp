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

/*
 * This program generated the files in this folder, it is kept here for reference (it is not built
 * by the Chaste build system).
 *
 * It was built against the ParameterBox.hpp/.cpp from ApPredict before bisection refinement was
 * introduced (git commit 0813da9), with Boost 1.74 (the oldest Boost used by ApPredict's CI and docker
 * images - text archives made with newer Boost versions can't be read by older ones).
 *
 * It creates a ParameterBox<2> refined by subdivision into 2^DIM boxes, as the LookupTableGenerator
 * used to do, archives it, and records interpolated values so that
 * TestLookupTableBackwardsCompatibility can check that newer code still loads and interpolates
 * these archives in exactly the same way.
 */
#include <fstream>
#include <iomanip>
#include "CheckpointArchiveTypes.hpp"
#include "Exception.hpp"
#include "ParameterBox.hpp"

double Separable2d(const c_vector<double, 2u>& rX)
{
    return exp(rX[0]) * (1.0 + rX[1] * rX[1]);
}

void AssignDataAtPoint(ParameterBox<2>& rBox, c_vector<double, 2u>* pPoint)
{
    std::vector<double> qoi(1u, Separable2d(*pPoint));
    rBox.AssignQoIValues(pPoint, boost::shared_ptr<ParameterPointData>(new ParameterPointData(qoi, 0u)));
}

int main()
{
    ParameterBox<2>* p_box = new ParameterBox<2>(NULL);
    std::set<c_vector<double, 2u>*, c_vector_compare<2u> > corners = p_box->GetCorners();
    for (auto iter = corners.begin(); iter != corners.end(); ++iter)
    {
        AssignDataAtPoint(*p_box, *iter);
    }

    // Refine as the old LookupTableGenerator did.
    const double tolerance = 5e-2;
    while (ParameterBox<2>* p_refine = p_box->FindBoxWithLargestQoIErrorEstimate(0u, tolerance))
    {
        std::set<c_vector<double, 2u>*, c_vector_compare<2u> > new_points = p_refine->SubDivide();
        for (auto iter = new_points.begin(); iter != new_points.end(); ++iter)
        {
            AssignDataAtPoint(*p_box, *iter);
        }
    }
    std::cout << "Created ParameterBox<2> with " << p_box->GetCorners().size() << " points.\n";

    {
        std::ofstream ofs("ParameterBox2d_legacy.arch");
        boost::archive::text_oarchive output_arch(ofs);
        output_arch << p_box;
    }

    std::ofstream interp_file("ParameterBox2d_legacy_interpolated.dat");
    interp_file << std::setprecision(17);
    for (unsigned i = 0; i <= 20u; i++)
    {
        for (unsigned j = 0; j <= 20u; j++)
        {
            c_vector<double, 2u> point;
            point[0] = 0.05 * i;
            point[1] = 0.03 + 0.047 * j;
            interp_file << point[0] << "\t" << point[1] << "\t" << p_box->InterpolateQoIsAt(point)[0] << "\n";
        }
    }
    delete p_box;
    return 0;
}
