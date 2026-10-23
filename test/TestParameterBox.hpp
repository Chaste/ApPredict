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

#ifndef TESTPARAMETERBOX_HPP_
#define TESTPARAMETERBOX_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

#include "OutputFileHandler.hpp"
#include "ParameterBox.hpp"

/** Helper to define a pointer to a test function without it being used for template argument deduction. */
template <unsigned DIM>
struct TestFunction
{
    typedef double (*Type)(const c_vector<double, DIM>&);
};

/**
 * Test the Parameter Box class.
 *
 */
class TestParameterBox : public CxxTest::TestSuite
{
private:
    /* Some cheap test functions to emulate QoIs with */
    static double Exp1d(const c_vector<double, 1u>& rX)
    {
        return exp(rX[0]);
    }

    static double ExpX2d(const c_vector<double, 2u>& rX)
    {
        return exp(rX[0]);
    }

    static double ExpTwoX2d(const c_vector<double, 2u>& rX)
    {
        return exp(2.0 * rX[0]);
    }

    static double Separable2d(const c_vector<double, 2u>& rX)
    {
        return exp(rX[0]) * (1.0 + rX[1] * rX[1]);
    }

    static double Bilinear2d(const c_vector<double, 2u>& rX)
    {
        return 1.0 + 2.0 * rX[0] - 3.0 * rX[1] + 4.0 * rX[0] * rX[1];
    }

    static double Trilinear3d(const c_vector<double, 3u>& rX)
    {
        return rX[0] + 2.0 * rX[1] + 3.0 * rX[2] + rX[0] * rX[1] * rX[2];
    }

    c_vector<double, 2u> MakePoint(double x, double y)
    {
        c_vector<double, 2u> point;
        point[0] = x;
        point[1] = y;
        return point;
    }

    // Assign data from a test function at a single point
    template <unsigned DIM>
    void AssignFunctionDataAtPoint(ParameterBox<DIM>& rBox,
                                   c_vector<double, DIM>* pPoint,
                                   typename TestFunction<DIM>::Type pFunction)
    {
        std::vector<double> qoi(1u, pFunction(*pPoint));
        boost::shared_ptr<ParameterPointData> p_data = boost::shared_ptr<ParameterPointData>(new ParameterPointData(qoi, 0u));
        rBox.AssignQoIValues(pPoint, p_data);
    }

    // Assign data from a test function at all the corners in a box
    template <unsigned DIM>
    void AssignFunctionData(ParameterBox<DIM>& rBox,
                            typename TestFunction<DIM>::Type pFunction)
    {
        std::vector<c_vector<double, DIM>*> corners = rBox.GetCornersAsVector();
        for (unsigned i = 0; i < corners.size(); i++)
        {
            AssignFunctionDataAtPoint<DIM>(rBox, corners[i], pFunction);
        }
    }

    // Mimic what the LookupTableGenerator does, but with a cheap function, returns number of points.
    template <unsigned DIM>
    unsigned RefineToTolerance(ParameterBox<DIM>& rBox,
                               typename TestFunction<DIM>::Type pFunction,
                               double tolerance,
                               bool bisect)
    {
        AssignFunctionData(rBox, pFunction);
        for (unsigned iteration = 0; iteration < 10000u; iteration++)
        {
            ParameterBox<DIM>* p_box = rBox.FindBoxWithLargestQoIErrorEstimate(0u, tolerance);
            if (!p_box)
            {
                break;
            }
            std::set<c_vector<double, DIM>*, c_vector_compare<DIM> > new_points;
            if (bisect)
            {
                new_points = p_box->SubDivide(p_box->ChooseDimensionToSplit(0u));
            }
            else
            {
                new_points = p_box->SubDivide();
            }
            for (auto iter = new_points.begin(); iter != new_points.end(); ++iter)
            {
                AssignFunctionDataAtPoint<DIM>(rBox, *iter, pFunction);
            }
        }
        return rBox.GetCornersAsVector().size();
    }

    template <unsigned DIM>
    std::vector<ParameterBox<DIM>*> GetLeafBoxes(ParameterBox<DIM>& rBox)
    {
        std::vector<ParameterBox<DIM>*> all_boxes = rBox.GetWholeFamilyOfBoxes();
        std::vector<ParameterBox<DIM>*> leaves;
        for (unsigned i = 0; i < all_boxes.size(); i++)
        {
            if (!all_boxes[i]->IsParent())
            {
                leaves.push_back(all_boxes[i]);
            }
        }
        return leaves;
    }

    // Largest interpolation error on a regular grid over [0,1]^2
    double MaxInterpolationError2d(ParameterBox<2>& rBox,
                                   double (*pFunction)(const c_vector<double, 2u>&))
    {
        double max_error = 0.0;
        for (unsigned i = 0; i <= 20u; i++)
        {
            for (unsigned j = 0; j <= 20u; j++)
            {
                c_vector<double, 2u> point = MakePoint(0.05 * i, 0.05 * j);
                double error = fabs(rBox.InterpolateQoIsAt(point)[0] - pFunction(point));
                max_error = std::max(max_error, error);
            }
        }
        return max_error;
    }

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

    void TestBisection1dMatchesFullSubdivision()
    {
        // In 1D a bisection is the same as a subdivision into 2^DIM boxes, so should give identical results.
        ParameterBox<1> legacy_box(NULL);
        ParameterBox<1> bisected_box(NULL);

        unsigned num_legacy_points = RefineToTolerance(legacy_box, Exp1d, 1e-4, false);
        unsigned num_bisected_points = RefineToTolerance(bisected_box, Exp1d, 1e-4, true);

        TS_ASSERT_EQUALS(num_legacy_points, num_bisected_points);
        TS_ASSERT_LESS_THAN(10u, num_bisected_points); // Check it did something

        std::vector<c_vector<double, 1u>*> legacy_points = legacy_box.GetCornersAsVector();
        std::vector<c_vector<double, 1u>*> bisected_points = bisected_box.GetCornersAsVector();
        for (unsigned i = 0; i < std::min(legacy_points.size(), bisected_points.size()); i++)
        {
            TS_ASSERT_DELTA((*legacy_points[i])[0], (*bisected_points[i])[0], 1e-15);
        }

        for (unsigned i = 0; i <= 100u; i++)
        {
            c_vector<double, 1u> point;
            point[0] = 0.01 * i;
            TS_ASSERT_DELTA(legacy_box.InterpolateQoIsAt(point)[0], bisected_box.InterpolateQoIsAt(point)[0], 1e-15);
        }
    }

    void TestBisection2dGeometryAndErrorEstimates()
    {
        ParameterBox<2> parent_box(NULL);
        AssignFunctionData(parent_box, ExpX2d);
        TS_ASSERT_EQUALS(parent_box.GetRefinementLevel(), 0u);

        // All dimensions have unknown errors, so tie-break on variation, which is only along x.
        TS_ASSERT_EQUALS(parent_box.ChooseDimensionToSplit(0u), 0u);

        TS_ASSERT_THROWS_THIS(parent_box.SubDivide(2u),
                              "Cannot subdivide along dimension 2 of a 2D box.");

        std::set<c_vector<double, 2u>*, c_vector_compare<2u> > new_points = parent_box.SubDivide(0u);
        TS_ASSERT_EQUALS(parent_box.IsParent(), true);
        TS_ASSERT_EQUALS(new_points.size(), 2u); // In 2D a bisection requires the addition of 2 new points.
        std::vector<c_vector<double, 2u>*> new_points_vec(new_points.begin(), new_points.end());
        TS_ASSERT_DELTA((*new_points_vec[0])[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*new_points_vec[0])[1], 0.0, 1e-12);
        TS_ASSERT_DELTA((*new_points_vec[1])[0], 0.5, 1e-12);
        TS_ASSERT_DELTA((*new_points_vec[1])[1], 1.0, 1e-12);

        TS_ASSERT_THROWS_THIS(parent_box.SubDivide(1u), "Already subdivided this box.");
        TS_ASSERT_THROWS_THIS(parent_box.ChooseDimensionToSplit(0u), "This box has already been subdivided.");
        TS_ASSERT_EQUALS(parent_box.GetCornersAsVector().size(), 6u);

        std::vector<ParameterBox<2>*> daughter_boxes = parent_box.GetDaughterBoxes();
        TS_ASSERT_EQUALS(daughter_boxes.size(), 2u);
        TS_ASSERT_EQUALS(daughter_boxes[0]->GetNewCorners().size(), 2u);
        TS_ASSERT_EQUALS(daughter_boxes[1]->GetNewCorners().size(), 0u);
        TS_ASSERT_DELTA(daughter_boxes[0]->mMin[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(daughter_boxes[0]->mMax[0], 0.5, 1e-12);
        TS_ASSERT_DELTA(daughter_boxes[0]->mMax[1], 1.0, 1e-12);
        TS_ASSERT_DELTA(daughter_boxes[1]->mMin[0], 0.5, 1e-12);
        TS_ASSERT_DELTA(daughter_boxes[1]->mMin[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(daughter_boxes[1]->mMax[0], 1.0, 1e-12);
        for (unsigned i = 0; i < 2u; i++)
        {
            TS_ASSERT_EQUALS(daughter_boxes[i]->GetRefinementLevel(), 1u);
            TS_ASSERT_EQUALS(daughter_boxes[i]->GetGeneration(), 1u);
            TS_ASSERT_EQUALS(daughter_boxes[i]->GetOwnCorners().size(), 4u);
            TS_ASSERT_THROWS_THIS(daughter_boxes[i]->GetMaxErrorsInPredictedQoIs(),
                                  "Not all the parameter points (which you can get with GetNewCorners()) have been assigned data. Error estimates unavailable.");
        }

        AssignFunctionData(parent_box, ExpX2d);

        // Linear interpolation between exp(0) and exp(1), compared to exp(0.5).
        const double expected_error = (exp(0.0) + exp(1.0)) / 2.0 - exp(0.5);
        for (unsigned i = 0; i < 2u; i++)
        {
            TS_ASSERT_DELTA(daughter_boxes[i]->GetMaxErrorsInPredictedQoIs()[0], expected_error, 1e-12);

            // We only know about errors in x, the error in y is unknown.
            std::vector<double> errors = daughter_boxes[i]->GetErrorEstimatesPerDimension(0u);
            TS_ASSERT_DELTA(errors[0], expected_error, 1e-12);
            TS_ASSERT_EQUALS(errors[1], DBL_MAX);
            TS_ASSERT_EQUALS(daughter_boxes[i]->GetMaxErrorInQoIEstimateInThisBox(0u), DBL_MAX);
            TS_ASSERT_EQUALS(daughter_boxes[i]->ChooseDimensionToSplit(0u), 1u);
        }

        // Split along y
        new_points = daughter_boxes[0]->SubDivide(1u);
        TS_ASSERT_EQUALS(new_points.size(), 2u); // (0,0.5) and (0.5,0.5)
        AssignFunctionData(parent_box, ExpX2d);

        ParameterBox<2>* p_grand_daughter = daughter_boxes[0]->GetDaughterBoxes()[0];
        TS_ASSERT_EQUALS(p_grand_daughter->GetRefinementLevel(), 2u);
        TS_ASSERT_DELTA(p_grand_daughter->mMax[0], 0.5, 1e-12);
        TS_ASSERT_DELTA(p_grand_daughter->mMax[1], 0.5, 1e-12);

        // exp(x) is constant in y so there is no error in the y direction,
        // and the x error estimate is inherited from the parent.
        TS_ASSERT_DELTA(p_grand_daughter->GetMaxErrorsInPredictedQoIs()[0], 0.0, 1e-12);
        std::vector<double> errors = p_grand_daughter->GetErrorEstimatesPerDimension(0u);
        TS_ASSERT_DELTA(errors[0], expected_error, 1e-12);
        TS_ASSERT_DELTA(errors[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_grand_daughter->GetMaxErrorInQoIEstimateInThisBox(0u), expected_error, 1e-12);
        TS_ASSERT_EQUALS(p_grand_daughter->ChooseDimensionToSplit(0u), 0u);

        // For comparison, a subdivision into 2^DIM boxes counts as DIM refinement levels,
        // and gives error estimates in every dimension.
        ParameterBox<2> legacy_box(NULL);
        AssignFunctionData(legacy_box, ExpX2d);
        legacy_box.SubDivide();
        AssignFunctionData(legacy_box, ExpX2d);
        ParameterBox<2>* p_legacy_daughter = legacy_box.GetDaughterBoxes()[0];
        TS_ASSERT_EQUALS(p_legacy_daughter->GetRefinementLevel(), 2u);
        errors = p_legacy_daughter->GetErrorEstimatesPerDimension(0u);
        TS_ASSERT_DELTA(errors[0], p_legacy_daughter->GetMaxErrorsInPredictedQoIs()[0], 1e-12);
        TS_ASSERT_DELTA(errors[1], p_legacy_daughter->GetMaxErrorsInPredictedQoIs()[0], 1e-12);
    }

    void TestBisectionWhenAllNewPointsAlreadyExist()
    {
        /*
         * Build up a situation where neighbouring boxes have already evaluated all the points
         * on the plane that a bisection needs, so no new points are required.
         */
        ParameterBox<2> parent_box(NULL);
        AssignFunctionData(parent_box, Separable2d);

        parent_box.SubDivide(0u); // [0,0.5]x[0,1] and [0.5,1]x[0,1]
        AssignFunctionData(parent_box, Separable2d);
        ParameterBox<2>* p_left = parent_box.GetDaughterBoxes()[0];

        p_left->SubDivide(1u); // [0,0.5]x[0,0.5] and [0,0.5]x[0.5,1]
        AssignFunctionData(parent_box, Separable2d);
        ParameterBox<2>* p_left_bottom = p_left->GetDaughterBoxes()[0];
        ParameterBox<2>* p_left_top = p_left->GetDaughterBoxes()[1];

        p_left_bottom->SubDivide(1u); // [0,0.5]x[0,0.25] and [0,0.5]x[0.25,0.5]
        AssignFunctionData(parent_box, Separable2d);
        ParameterBox<2>* p_lbb = p_left_bottom->GetDaughterBoxes()[0];
        ParameterBox<2>* p_lbt = p_left_bottom->GetDaughterBoxes()[1];

        TS_ASSERT_EQUALS(p_lbb->SubDivide(0u).size(), 2u); // Creates (0.25,0) and (0.25,0.25)
        AssignFunctionData(parent_box, Separable2d);
        TS_ASSERT_EQUALS(p_left_top->SubDivide(0u).size(), 2u); // Creates (0.25,0.5) and (0.25,1)
        AssignFunctionData(parent_box, Separable2d);

        // Now both (0.25,0.25) and (0.25,0.5) have been evaluated already.
        std::set<c_vector<double, 2u>*, c_vector_compare<2u> > new_points = p_lbt->SubDivide(0u);
        TS_ASSERT_EQUALS(new_points.size(), 0u);

        // But we should still have error estimates, from interpolating across p_lbt.
        double error_low = fabs(0.5 * (Separable2d(MakePoint(0.0, 0.25)) + Separable2d(MakePoint(0.5, 0.25)))
                                - Separable2d(MakePoint(0.25, 0.25)));
        double error_high = fabs(0.5 * (Separable2d(MakePoint(0.0, 0.5)) + Separable2d(MakePoint(0.5, 0.5)))
                                 - Separable2d(MakePoint(0.25, 0.5)));
        std::vector<ParameterBox<2>*> daughters = p_lbt->GetDaughterBoxes();
        TS_ASSERT_EQUALS(daughters.size(), 2u);
        for (unsigned i = 0; i < daughters.size(); i++)
        {
            TS_ASSERT_EQUALS(daughters[i]->mAllCornersEvaluated, true);
            TS_ASSERT_EQUALS(daughters[i]->GetMaxErrorsInPredictedQoIs().size(), 1u);
            TS_ASSERT_DELTA(daughters[i]->GetMaxErrorsInPredictedQoIs()[0], std::max(error_low, error_high), 1e-12);
            TS_ASSERT_DELTA(daughters[i]->GetErrorEstimatesPerDimension(0u)[0], std::max(error_low, error_high), 1e-12);
        }

        // And the rest of the tree should still work as usual.
        TS_ASSERT(parent_box.FindBoxWithLargestQoIErrorEstimate(0u, 1e-6) != NULL);
        double percentage = parent_box.ReportPercentageOfSpaceWhereToleranceIsMetForQoI(1.0, 0u);
        TS_ASSERT_LESS_THAN_EQUALS(0.0, percentage);
        TS_ASSERT_LESS_THAN_EQUALS(percentage, 100.0);

        // Interpolation in the refined region [0,0.5)x[0.25,0.5] should be reasonable.
        // (Points on x=0.5 are interpolated in the unrefined box to the right.)
        for (unsigned i = 0; i < 10u; i++)
        {
            for (unsigned j = 0; j <= 10u; j++)
            {
                c_vector<double, 2u> point = MakePoint(0.05 * i, 0.25 + 0.025 * j);
                TS_ASSERT_DELTA(parent_box.InterpolateQoIsAt(point)[0], Separable2d(point), 5e-2);
            }
        }
    }

    void TestBisection3dVisitsEveryDimensionOnce()
    {
        // A multilinear function has no interpolation error, but we need to split
        // each part of the tree along each dimension once to find that out.
        // So we should end up with the same 3x3x3 grid as one subdivision into 2^DIM.
        ParameterBox<3> parent_box(NULL);
        TS_ASSERT_EQUALS(RefineToTolerance(parent_box, Trilinear3d, 1e-6, true), 27u);

        std::vector<ParameterBox<3>*> leaves = GetLeafBoxes(parent_box);
        TS_ASSERT_EQUALS(leaves.size(), 8u);
        for (unsigned i = 0; i < leaves.size(); i++)
        {
            TS_ASSERT_EQUALS(leaves[i]->GetRefinementLevel(), 3u);
            for (unsigned j = 0; j < 3u; j++)
            {
                TS_ASSERT_DELTA(leaves[i]->mMax[j] - leaves[i]->mMin[j], 0.5, 1e-12);
            }
        }

        for (unsigned i = 0; i <= 10u; i++)
        {
            c_vector<double, 3u> point;
            point[0] = 0.1 * i;
            point[1] = 0.07 * i;
            point[2] = 1.0 - 0.1 * i;
            TS_ASSERT_DELTA(parent_box.InterpolateQoIsAt(point)[0], Trilinear3d(point), 1e-12);
        }
    }

    void TestBisectionOnlyRefinesDimensionsThatNeedIt()
    {
        const double tolerance = 1e-2;

        // exp(2x) only varies in x
        ParameterBox<2> bisected_box(NULL);
        ParameterBox<2> legacy_box(NULL);
        unsigned num_bisected_points = RefineToTolerance(bisected_box, ExpTwoX2d, tolerance, true);
        unsigned num_legacy_points = RefineToTolerance(legacy_box, ExpTwoX2d, tolerance, false);
        std::cout << "exp(2x) to tolerance " << tolerance << ": bisection used " << num_bisected_points
                  << " points, subdivision into 2^DIM used " << num_legacy_points << " points.\n";
        TS_ASSERT_LESS_THAN(num_bisected_points, num_legacy_points);

        // We only needed to split once in y, all refinement after that is in x.
        std::vector<ParameterBox<2>*> leaves = GetLeafBoxes(bisected_box);
        for (unsigned i = 0; i < leaves.size(); i++)
        {
            TS_ASSERT_DELTA(leaves[i]->mMax[1] - leaves[i]->mMin[1], 0.5, 1e-12);
        }

        TS_ASSERT_LESS_THAN(MaxInterpolationError2d(bisected_box, ExpTwoX2d), tolerance);
        TS_ASSERT_LESS_THAN(MaxInterpolationError2d(legacy_box, ExpTwoX2d), tolerance);
        TS_ASSERT_DELTA(bisected_box.ReportPercentageOfSpaceWhereToleranceIsMetForQoI(tolerance, 0u), 100.0, 1e-9);

        // A function that is curved in both dimensions
        ParameterBox<2> bisected_box2(NULL);
        ParameterBox<2> legacy_box2(NULL);
        num_bisected_points = RefineToTolerance(bisected_box2, Separable2d, tolerance, true);
        num_legacy_points = RefineToTolerance(legacy_box2, Separable2d, tolerance, false);
        std::cout << "exp(x)(1+y^2) to tolerance " << tolerance << ": bisection used " << num_bisected_points
                  << " points, subdivision into 2^DIM used " << num_legacy_points << " points.\n";
        TS_ASSERT_LESS_THAN(num_bisected_points, num_legacy_points);
        TS_ASSERT_LESS_THAN(MaxInterpolationError2d(bisected_box2, Separable2d), tolerance);
    }

    void TestInterpolationOnMixedTree()
    {
        // A bilinear function should be interpolated exactly on any tree of rectangles,
        // here we mix up subdivisions into 2^DIM with bisections.
        ParameterBox<2> parent_box(NULL);
        AssignFunctionData(parent_box, Bilinear2d);
        parent_box.SubDivide();
        AssignFunctionData(parent_box, Bilinear2d);
        std::vector<ParameterBox<2>*> daughters = parent_box.GetDaughterBoxes();
        TS_ASSERT_EQUALS(daughters.size(), 4u);

        daughters[0]->SubDivide(0u);
        AssignFunctionData(parent_box, Bilinear2d);
        daughters[3]->SubDivide(1u);
        AssignFunctionData(parent_box, Bilinear2d);
        daughters[3]->GetDaughterBoxes()[1]->SubDivide(0u);
        AssignFunctionData(parent_box, Bilinear2d);
        daughters[1]->SubDivide();
        AssignFunctionData(parent_box, Bilinear2d);

        TS_ASSERT_DELTA(MaxInterpolationError2d(parent_box, Bilinear2d), 0.0, 1e-12);

        // Every leaf has had every dimension split somewhere in its family tree,
        // and there is no error, so nothing needs refining.
        std::vector<ParameterBox<2>*> leaves = GetLeafBoxes(parent_box);
        TS_ASSERT_EQUALS(leaves.size(), 10u);
        for (unsigned i = 0; i < leaves.size(); i++)
        {
            std::vector<double> errors = leaves[i]->GetErrorEstimatesPerDimension(0u);
            TS_ASSERT_DELTA(errors[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(errors[1], 0.0, 1e-12);
        }
        TS_ASSERT(parent_box.FindBoxWithLargestQoIErrorEstimate(0u, 1e-10) == NULL);
        TS_ASSERT_DELTA(parent_box.ReportPercentageOfSpaceWhereToleranceIsMetForQoI(1e-10, 0u), 100.0, 1e-9);
    }

    void TestBisectingATreeMadeBySubdivision()
    {
        // This is what happens when an old lookup table is refined further.
        ParameterBox<2> parent_box(NULL);
        unsigned num_legacy_points = RefineToTolerance(parent_box, Separable2d, 5e-2, false);
        TS_ASSERT_LESS_THAN(MaxInterpolationError2d(parent_box, Separable2d), 5e-2);

        unsigned num_points = RefineToTolerance(parent_box, Separable2d, 5e-3, true);
        TS_ASSERT_LESS_THAN(num_legacy_points, num_points);
        TS_ASSERT_LESS_THAN(MaxInterpolationError2d(parent_box, Separable2d), 5e-3);
        TS_ASSERT_DELTA(parent_box.ReportPercentageOfSpaceWhereToleranceIsMetForQoI(5e-3, 0u), 100.0, 1e-9);
    }

    void TestArchivingBisectedParameterBox()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "BisectedParameterBox.arch";

        const double tolerance = 1e-2;
        unsigned num_points;
        std::vector<double> interpolated_values;
        c_vector<double, 2u> next_box_min;
        unsigned next_dimension;

        // SAVE
        {
            ParameterBox<2>* p_box = new ParameterBox<2>(NULL);
            num_points = RefineToTolerance(*p_box, Separable2d, tolerance, true);
            for (unsigned i = 0; i <= 20u; i++)
            {
                interpolated_values.push_back(p_box->InterpolateQoIsAt(MakePoint(0.05 * i, 1.0 - 0.04 * i))[0]);
            }
            ParameterBox<2>* p_next_box = p_box->FindBoxWithLargestQoIErrorEstimate(0u, 0.5 * tolerance);
            TS_ASSERT(p_next_box);
            next_box_min = p_next_box->mMin;
            next_dimension = p_next_box->ChooseDimensionToSplit(0u);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_box;
            delete p_box;
        }

        // LOAD
        {
            ParameterBox<2>* p_box;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_box;

            TS_ASSERT_EQUALS(p_box->GetCornersAsVector().size(), num_points);
            for (unsigned i = 0; i <= 20u; i++)
            {
                TS_ASSERT_DELTA(p_box->InterpolateQoIsAt(MakePoint(0.05 * i, 1.0 - 0.04 * i))[0], interpolated_values[i], 1e-12);
            }
            ParameterBox<2>* p_next_box = p_box->FindBoxWithLargestQoIErrorEstimate(0u, 0.5 * tolerance);
            TS_ASSERT(p_next_box);
            TS_ASSERT_DELTA(p_next_box->mMin[0], next_box_min[0], 1e-12);
            TS_ASSERT_DELTA(p_next_box->mMin[1], next_box_min[1], 1e-12);
            TS_ASSERT_EQUALS(p_next_box->ChooseDimensionToSplit(0u), next_dimension);

            // And we can carry on refining
            TS_ASSERT_LESS_THAN(num_points, RefineToTolerance(*p_box, Separable2d, 0.5 * tolerance, true));

            delete p_box;
        }
    }
};

#endif // TESTPARAMETERBOX_HPP_
