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

#ifndef TESTPARAMETERBOX_HPP_
#define TESTPARAMETERBOX_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

#include "OutputFileHandler.hpp"
#include "ParameterBox.hpp"

/**
 * Test the Parameter Box class.
 *
 */
class TestParameterBox : public CxxTest::TestSuite
{
private:
    // 1D box
    void AssignExponentialData(ParameterBox<1>& rBox,
                               std::vector<c_vector<double, 1u>*>& rCorners)
    {
        unsigned error_code = 0u;
        // Invent some initial guesses (some too big, some too small)
        for (unsigned i = 0; i < rCorners.size(); i++)
        {
            std::vector<double> qoi;
            qoi.push_back(0.5); // Our initial guess is 0.5 everywhere.

            boost::shared_ptr<ParameterPointData> p_data = boost::shared_ptr<ParameterPointData>(new ParameterPointData(qoi, error_code));

            // Assign this data as an estimate at this corner (with the 'true' flag).
            rBox.AssignQoIValues(rCorners[i], p_data, true);

            // Check that the error estimates are being assigned to the parameter point data appropriately.
            TS_ASSERT_EQUALS(p_data->HasErrorEstimates(), false);
            TS_ASSERT_THROWS_THIS(p_data->rGetQoIErrorEstimates(),
                                  "Error estimates have not been set on this parameter data point.");
        }

        // Assign some 'real' data with exp(x)
        for (unsigned i = 0; i < rCorners.size(); i++)
        {
            std::vector<double> qoi;
            qoi.push_back(exp((*(rCorners[i]))[0]));
            boost::shared_ptr<ParameterPointData> p_data = boost::shared_ptr<ParameterPointData>(new ParameterPointData(qoi, error_code));
            rBox.AssignQoIValues(rCorners[i], p_data);
        }
    }

    // 2D Box just the same (exponential data on x only)
    void AssignExponentialData(ParameterBox<2>& rBox,
                               std::vector<c_vector<double, 2u>*>& rCorners)
    {
        unsigned error_code = 0u;
        // Invent some data with exp(x)
        for (unsigned i = 0; i < rCorners.size(); i++)
        {
            std::vector<double> qoi;
            qoi.push_back(exp((*(rCorners[i]))[0]));
            boost::shared_ptr<ParameterPointData> p_data = boost::shared_ptr<ParameterPointData>(new ParameterPointData(qoi, error_code));
            rBox.AssignQoIValues(rCorners[i], p_data);
        }
    }

public:
    void TestParameterBox1d()
    {
        ParameterBox<1> parent_box_1d(NULL);

        TS_ASSERT_EQUALS(parent_box_1d.GetGeneration(), 0u);
        std::vector<c_vector<double, 1u>*> corner_parameters = parent_box_1d.GetCornersAsVector();

        TS_ASSERT_EQUALS(corner_parameters.size(), 2u);
        TS_ASSERT_DELTA((*(corner_parameters[0]))[0], 0.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[1]))[0], 1.0, 1e-12);

        TS_ASSERT_THROWS_THIS(parent_box_1d.GetMaxErrorsInPredictedQoIs(),
                              "Not all the parameter points (which you can get with GetNewCorners()) have been assigned data. Error estimates unavailable.");

        AssignExponentialData(parent_box_1d, corner_parameters);

        c_vector<double, 1u> location;
        location[0] = 1.1;
        TS_ASSERT_THROWS_THIS(parent_box_1d.GetBoxContainingPoint(location),
                              "This point is not contained within this box (or any of its children).");

        location[0] = 0.44;
        ParameterBox<1>* p_box = parent_box_1d.GetBoxContainingPoint(location);
        TS_ASSERT_EQUALS(p_box, &parent_box_1d);

        // Check some exceptions
        TS_ASSERT_THROWS_THIS(parent_box_1d.GetParent(),
                              "This parameter box has no parent.");

        // Check division of the box
        TS_ASSERT_EQUALS(parent_box_1d.IsParent(), false);
        std::set<c_vector<double, 1u>*, c_vector_compare<1u> > new_points = parent_box_1d.SubDivide();
        TS_ASSERT_EQUALS(parent_box_1d.IsParent(), true);
        TS_ASSERT_EQUALS(new_points.size(), 1u); // In 1D a SubDivide requires the addition of 1 new point.

        TS_ASSERT_THROWS_THIS(parent_box_1d.SubDivide(),
                              "Already subdivided this box.");

        std::vector<ParameterBox<1>*> daughter_boxes = parent_box_1d.GetDaughterBoxes();
        TS_ASSERT_EQUALS(daughter_boxes.size(), 2u);
        TS_ASSERT_EQUALS(daughter_boxes[0]->GetGeneration(), 1u);
        TS_ASSERT_EQUALS(daughter_boxes[1]->GetGeneration(), 1u);

        // Box 0
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[0]))[0], 0.0, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[1]))[0], 0.5, 1e-12);

        std::set<c_vector<double, 1u>*, c_vector_compare<1u> > new_corners = daughter_boxes[0]->GetNewCorners();
        TS_ASSERT_EQUALS(new_corners.size(), 1u);
        double new_corner_pos = (**(new_corners.begin()))[0];
        TS_ASSERT_DELTA(new_corner_pos, 0.5, 1e-12);

        // Box 1
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[0]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[1]))[0], 1.0, 1e-12);

        std::set<c_vector<double, 1u>*, c_vector_compare<1u> > new_corners2 = daughter_boxes[1]->GetNewCorners();
        TS_ASSERT_EQUALS(new_corners2.size(), 0u);

        corner_parameters = parent_box_1d.GetCornersAsVector();
        TS_ASSERT_EQUALS(corner_parameters.size(), 3u);
        TS_ASSERT_DELTA((*(corner_parameters[0]))[0], 0.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[1]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[2]))[0], 1.0, 1e-12);

        p_box = parent_box_1d.GetBoxContainingPoint(location);
        TS_ASSERT_EQUALS(p_box, daughter_boxes[0]);

        // Check nesting works
        AssignExponentialData(parent_box_1d, corner_parameters);

        std::vector<double> errors1 = daughter_boxes[0]->GetMaxErrorsInPredictedQoIs();
        std::vector<double> errors2 = daughter_boxes[1]->GetMaxErrorsInPredictedQoIs();

        /////////////////////////////////////////////////////
        // TEST THE INTEPOLATION AND ERROR ESTIMATION SCHEME
        ////////////////////////////////////////////////////

        // Since these boxes share the same new point, this should be the same number
        TS_ASSERT_DELTA(errors1[0], errors2[0], 1e-12);

        // The error should be linear approximation between exp(0) and exp(1),
        // compared to exp(0.5). This is very pleasing!
        TS_ASSERT_DELTA(errors1[0], (exp(0) + exp(1)) / 2.0 - exp(0.5), 1e-12);

        // Divide box at 0,0.5 into two:
        daughter_boxes[0]->SubDivide();

        corner_parameters = parent_box_1d.GetCornersAsVector();
        TS_ASSERT_EQUALS(corner_parameters.size(), 4u);
        TS_ASSERT_DELTA((*(corner_parameters[0]))[0], 0.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[1]))[0], 0.25, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[2]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[3]))[0], 1.0, 1e-12);

        // And the other box:
        AssignExponentialData(parent_box_1d, corner_parameters);
        daughter_boxes[1]->SubDivide();

        corner_parameters = parent_box_1d.GetCornersAsVector();
        TS_ASSERT_EQUALS(corner_parameters.size(), 5u);
        TS_ASSERT_DELTA((*(corner_parameters[0]))[0], 0.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[1]))[0], 0.25, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[2]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[3]))[0], 0.75, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[4]))[0], 1.0, 1e-12);

        AssignExponentialData(parent_box_1d, corner_parameters);

        // The largest error in exp(x) estimates should be at 0.75.
        // So either box 0.5->0.75 or 0.75->1.0 would be fine.
        p_box = parent_box_1d.FindBoxWithLargestQoIErrorEstimate(0u, DBL_MIN);

        TS_ASSERT(p_box);
        TS_ASSERT_EQUALS(p_box->IsParent(), false);
        TS_ASSERT_EQUALS(p_box->GetGeneration(), 2u);
        TS_ASSERT_DELTA(p_box->GetMaxErrorsInPredictedQoIs()[0],
                        (exp(1) - exp(0.5)) / 2.0 + exp(0.5) - exp(0.75), 1e-12);
        TS_ASSERT_DELTA((*(p_box->GetCornersAsVector()[0]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(p_box->GetCornersAsVector()[1]))[0], 0.75, 1e-12);

        // Test the self-interpolation methods.
        location[0] = 0.0;
        std::vector<double> interp = parent_box_1d.InterpolateQoIsAt(location);
        TS_ASSERT_EQUALS(interp.size(), 1.0);
        TS_ASSERT_DELTA(interp[0], 1.0, 1e-12);

        location[0] = 1.0;
        interp = parent_box_1d.InterpolateQoIsAt(location);
        TS_ASSERT_DELTA(interp[0], exp(1.0), 1e-12);

        location[0] = 0.44;
        interp = parent_box_1d.InterpolateQoIsAt(location);
        TS_ASSERT_DELTA(interp[0], exp(0.44), 1e-2);
    }

    void TestArchivingParameterBox()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ParameterBox.arch";

        // SAVE
        {
            // Repeat most of the first test to get us to a state with 5 boxes...

            ParameterBox<1>* p_parent_box_1d = new ParameterBox<1>(NULL);
            std::vector<c_vector<double, 1u>*> corner_parameters = p_parent_box_1d->GetCornersAsVector();

            TS_ASSERT_EQUALS(corner_parameters.size(), 2u);

            AssignExponentialData(*p_parent_box_1d, corner_parameters);
            p_parent_box_1d->SubDivide();

            std::vector<ParameterBox<1>*> daughter_boxes = p_parent_box_1d->GetDaughterBoxes();
            TS_ASSERT_EQUALS(daughter_boxes.size(), 2u);

            std::set<c_vector<double, 1u>*, c_vector_compare<1u> > new_corners = daughter_boxes[0]->GetNewCorners();
            TS_ASSERT_EQUALS(new_corners.size(), 1u);
            TS_ASSERT_DELTA((**(new_corners.begin()))[0], 0.5, 1e-12);

            new_corners = daughter_boxes[1]->GetNewCorners();
            TS_ASSERT_EQUALS(new_corners.size(), 0u);

            corner_parameters = p_parent_box_1d->GetCornersAsVector();
            TS_ASSERT_EQUALS(corner_parameters.size(), 3u);

            // Check nesting works
            // Divide box at 0,0.5 into two:
            AssignExponentialData(*p_parent_box_1d, corner_parameters);
            daughter_boxes[0]->SubDivide();

            corner_parameters = p_parent_box_1d->GetCornersAsVector();
            TS_ASSERT_EQUALS(corner_parameters.size(), 4u);

            // And the other box:
            AssignExponentialData(*p_parent_box_1d, corner_parameters);
            daughter_boxes[1]->SubDivide();

            corner_parameters = p_parent_box_1d->GetCornersAsVector();
            TS_ASSERT_EQUALS(corner_parameters.size(), 5u);

            AssignExponentialData(*p_parent_box_1d, corner_parameters);

            // The largest jump in exp(x) should be between 0.75 and 1.0
            ParameterBox<1u>* p_best_box = p_parent_box_1d->FindBoxWithLargestQoIErrorEstimate(0u, DBL_MIN);

            TS_ASSERT(p_best_box);
            TS_ASSERT_EQUALS(p_best_box->IsParent(), false);
            TS_ASSERT_DELTA(p_best_box->GetMaxErrorsInPredictedQoIs()[0],
                            (exp(1) - exp(0.5)) / 2.0 + exp(0.5) - exp(0.75), 1e-12);
            TS_ASSERT_DELTA((*(p_best_box->GetCornersAsVector()[0]))[0], 0.5, 1e-12);
            TS_ASSERT_DELTA((*(p_best_box->GetCornersAsVector()[1]))[0], 0.75, 1e-12);

            daughter_boxes = p_parent_box_1d->GetDaughterBoxes();
            TS_ASSERT_EQUALS(daughter_boxes.size(), 2u);

            std::vector<ParameterBox<1u>*> all_boxes = p_parent_box_1d->GetWholeFamilyOfBoxes();
            TS_ASSERT_EQUALS(all_boxes.size(), 7u);

            // Archive it

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_parent_box_1d;

            // Clean up memory
            delete p_parent_box_1d;
        }

        // LOAD
        {
            ParameterBox<1>* p_box;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> p_box;

            std::vector<ParameterBox<1>*> daughter_boxes = p_box->GetDaughterBoxes();
            TS_ASSERT_EQUALS(daughter_boxes.size(), 2u);

            ParameterBox<1u>* p_best_box = p_box->FindBoxWithLargestQoIErrorEstimate(0u, DBL_MIN);

            TS_ASSERT(p_best_box);
            TS_ASSERT_EQUALS(p_best_box->IsParent(), false);
            TS_ASSERT_DELTA(p_best_box->GetMaxErrorsInPredictedQoIs()[0],
                            (exp(1) - exp(0.5)) / 2.0 + exp(0.5) - exp(0.75), 1e-12);
            TS_ASSERT_DELTA((*(p_best_box->GetCornersAsVector()[0]))[0], 0.5, 1e-12);
            TS_ASSERT_DELTA((*(p_best_box->GetCornersAsVector()[1]))[0], 0.75, 1e-12);

            // Clean up memory, usually done by LookupTableGenerator
            delete p_box;
        }
    }

    void TestParameterBox2d()
    {
        ParameterBox<2> parent_box_2d(NULL);
        std::vector<c_vector<double, 2u>*> corner_parameters = parent_box_2d.GetCornersAsVector();

        TS_ASSERT_EQUALS(corner_parameters.size(), 4u);
        TS_ASSERT_DELTA((*(corner_parameters[0]))[0], 0.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[0]))[1], 0.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[1]))[0], 0.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[1]))[1], 1.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[2]))[0], 1.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[2]))[1], 0.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[3]))[0], 1.0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[3]))[1], 1.0, 1e-12);

        AssignExponentialData(parent_box_2d, corner_parameters);

        TS_ASSERT_EQUALS(parent_box_2d.IsParent(), false);
        std::set<c_vector<double, 2u>*, c_vector_compare<2u> > new_points = parent_box_2d.SubDivide();
        TS_ASSERT_EQUALS(parent_box_2d.IsParent(), true);
        TS_ASSERT_EQUALS(new_points.size(), 5u); // In 2D a SubDivide requires the addition of 5 new points.

        corner_parameters = parent_box_2d.GetCornersAsVector();
        TS_ASSERT_EQUALS(corner_parameters.size(), 9u);

        TS_ASSERT_DELTA((*(corner_parameters[0]))[0], 0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[0]))[1], 0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[1]))[0], 0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[1]))[1], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[2]))[0], 0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[2]))[1], 1, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[3]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[3]))[1], 0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[4]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[4]))[1], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[5]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[5]))[1], 1, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[6]))[0], 1, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[6]))[1], 0, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[7]))[0], 1, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[7]))[1], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[8]))[0], 1, 1e-12);
        TS_ASSERT_DELTA((*(corner_parameters[8]))[1], 1, 1e-12);

        std::vector<ParameterBox<2>*> daughter_boxes = parent_box_2d.GetDaughterBoxes();
        TS_ASSERT_EQUALS(daughter_boxes.size(), 4u);

        // Box 0
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[0]))[0], 0, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[0]))[1], 0, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[1]))[0], 0, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[1]))[1], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[2]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[2]))[1], 0, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[3]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[0]->GetCornersAsVector()[3]))[1], 0.5, 1e-12);
        // Box 1
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[0]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[0]))[1], 0, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[1]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[1]))[1], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[2]))[0], 1, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[2]))[1], 0, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[3]))[0], 1, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[1]->GetCornersAsVector()[3]))[1], 0.5, 1e-12);
        // Box 2
        TS_ASSERT_DELTA((*(daughter_boxes[2]->GetCornersAsVector()[0]))[0], 0, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[2]->GetCornersAsVector()[0]))[1], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[2]->GetCornersAsVector()[1]))[0], 0, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[2]->GetCornersAsVector()[1]))[1], 1, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[2]->GetCornersAsVector()[2]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[2]->GetCornersAsVector()[2]))[1], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[2]->GetCornersAsVector()[3]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[2]->GetCornersAsVector()[3]))[1], 1, 1e-12);
        // Box 3
        TS_ASSERT_DELTA((*(daughter_boxes[3]->GetCornersAsVector()[0]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[3]->GetCornersAsVector()[0]))[1], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[3]->GetCornersAsVector()[1]))[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[3]->GetCornersAsVector()[1]))[1], 1, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[3]->GetCornersAsVector()[2]))[0], 1, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[3]->GetCornersAsVector()[2]))[1], 0.5, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[3]->GetCornersAsVector()[3]))[0], 1, 1e-12);
        TS_ASSERT_DELTA((*(daughter_boxes[3]->GetCornersAsVector()[3]))[1], 1, 1e-12);
    }
};

#endif // TESTPARAMETERBOX_HPP_
