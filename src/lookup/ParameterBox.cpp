/*

Copyright (c) 2005-2017, University of Oxford.
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

#include <bitset> // for binary ops.

#include "Exception.hpp"
#include "ParameterBox.hpp"

template <unsigned DIM>
ParameterBox<DIM>::ParameterBox(ParameterBox<DIM>* pParent,
                                const c_vector<double, DIM>& rMin,
                                const c_vector<double, DIM>& rMax)
        : mAmParent(false),
          mpParentBox(pParent),
          mMin(rMin),
          mMax(rMax),
          mAllCornersEvaluated(false)
{
    if (!mpParentBox)
    {
        // I am the great grand parent box.
        mGeneration = 0u;
        mpGreatGrandParentBox = this;
    }
    else
    {
        mGeneration = mpParentBox->GetGeneration() + 1u;
        mpGreatGrandParentBox = mpParentBox->GetGreatGrandParent();
    }

    c_vector<double, DIM>* p_existing_corner; // Some working memory

    // Create a basic grid of the corners of this new box first
    for (unsigned i = 0; i < pow(2, DIM); i++)
    {
        // Use a binary conversion to get the right indices in place for the corners
        std::bitset<DIM> bin_i(i);
        c_vector<double, DIM> new_corner;
        for (unsigned j = 0; j < DIM; j++)
        {
            new_corner[j] = mMin[j] + (mMax[j] - mMin[j]) * (double)(bin_i[j]);
        }

        // Make a pointer to new location
        c_vector<double, DIM>* p_new_corner = new c_vector<double, DIM>(new_corner);

        // See if any of the other boxes have a corner that matches this location,
        // Because we made a smart iterator with our own comparison method this should work!
        DataMapIter iter = mpGreatGrandParentBox->mParameterPointDataMap.find(p_new_corner);
        if (iter != mpGreatGrandParentBox->mParameterPointDataMap.end())
        {
            // Also take note of any QoIs that have been evaluated on parent.
            p_existing_corner = (*iter).first;
            mCorners.push_back(p_existing_corner);
            mParameterPointDataMap[p_existing_corner] = (*iter).second;

            // If the existing corner is in a neighbouring box from this sub-division,
            // and therefore has not been evaluated yet, also prepare a prediction site here.
            if ((*iter).second == boost::shared_ptr<ParameterPointData>())
            {
                mParameterPointDataMapPredictions[p_existing_corner] = boost::shared_ptr<ParameterPointData>();
            }
            delete p_new_corner;
        }
        else
        {
            // Make a new corner map here
            mParameterPointDataMap[p_new_corner] = boost::shared_ptr<ParameterPointData>();
            mParameterPointDataMapPredictions[p_new_corner] = boost::shared_ptr<ParameterPointData>();
            // and store the real data points also on the great grandparent
            mpGreatGrandParentBox->mParameterPointDataMap[p_new_corner] = boost::shared_ptr<ParameterPointData>();
            // Make a list of new box corner
            mCorners.push_back(p_new_corner);
            mNewCorners.insert(p_new_corner);
        }
    }
}

template <unsigned DIM>
ParameterBox<DIM>::~ParameterBox()
{
    // If I am the original parent box, get all the corners,
    // and delete them all. Nobody else does anything.
    if (!mpParentBox)
    {
        for (DataMapIter iter = mParameterPointDataMap.begin();
             iter != mParameterPointDataMap.end();
             ++iter)
        {
            delete (*iter).first; // Delete the corner pointer
        }
    }

    // Also delete any daughter boxes that we made.
    if (mAmParent)
    {
        for (unsigned i = 0; i < mDaughterBoxes.size(); i++)
        {
            delete mDaughterBoxes[i];
        }
    }
}

template <unsigned DIM>
std::set<c_vector<double, DIM>*, c_vector_compare<DIM> > ParameterBox<DIM>::GetNewCorners()
{
    return mNewCorners;
}

template <unsigned DIM>
unsigned ParameterBox<DIM>::GetGeneration()
{
    return mGeneration;
}

template <unsigned DIM>
std::set<c_vector<double, DIM>*, c_vector_compare<DIM> > ParameterBox<DIM>::GetCorners()
{
    CornerSet corner_set;

    // Get any corners that we own,
    for (unsigned i = 0; i < mCorners.size(); i++)
    {
        corner_set.insert(mCorners[i]);
    }

    // Also include any corners that our children own.
    for (unsigned i = 0; i < mDaughterBoxes.size(); i++)
    {
        // For each parameter box
        CornerSet corners = mDaughterBoxes[i]->GetCorners();
        corner_set.insert(corners.begin(), corners.end());
    }

    return corner_set;
}

template <unsigned DIM>
std::vector<c_vector<double, DIM>*> ParameterBox<DIM>::GetCornersAsVector()
{
    CornerSet corner_set = GetCorners();
    // Convert to a std::vector
    std::vector<c_vector<double, DIM>*> result(corner_set.begin(), corner_set.end());
    return result;
}

template <unsigned DIM>
std::vector<c_vector<double, DIM>*> ParameterBox<DIM>::GetOwnCorners()
{
    return mCorners;
}

template <unsigned DIM>
bool ParameterBox<DIM>::IsParent()
{
    return mAmParent;
}

template <unsigned DIM>
ParameterBox<DIM>* ParameterBox<DIM>::GetParent()
{
    if (!mpParentBox)
    {
        EXCEPTION("This parameter box has no parent.");
    }
    return mpParentBox;
}

template <unsigned DIM>
ParameterBox<DIM>* ParameterBox<DIM>::GetGreatGrandParent()
{
    return mpGreatGrandParentBox;
}

template <unsigned DIM>
std::vector<ParameterBox<DIM>*> ParameterBox<DIM>::GetDaughterBoxes()
{
    return mDaughterBoxes;
}

template <unsigned DIM>
std::vector<ParameterBox<DIM>*> ParameterBox<DIM>::GetWholeFamilyOfBoxes()
{
    std::vector<ParameterBox<DIM>*> all_boxes;
    for (unsigned i = 0; i < mDaughterBoxes.size(); i++)
    {
        std::vector<ParameterBox<DIM>*> daughter_family_boxes = mDaughterBoxes[i]->GetWholeFamilyOfBoxes();
        for (unsigned j = 0; j < daughter_family_boxes.size(); j++)
        {
            all_boxes.push_back(daughter_family_boxes[j]);
        }
    }
    all_boxes.push_back(this);
    return all_boxes;
}

template <unsigned DIM>
std::set<c_vector<double, DIM>*, c_vector_compare<DIM> > ParameterBox<DIM>::SubDivide()
{
    if (mAmParent)
    {
        EXCEPTION("Already subdivided this box.");
    }

    CornerSet new_corners;

    // In order to subdivide a box we require there to be data in existence at all its corners.
    assert(mParameterPointDataMap.size() == mCorners.size());

    // Loop over each daughter box, we need to create them all at once.
    for (unsigned i = 0; i < pow(2, DIM); i++)
    {
        // Use a binary conversion to get the right indices in place.
        std::bitset<DIM> bin_i(i);
        c_vector<double, DIM> corner;
        for (unsigned j = 0; j < DIM; j++)
        {
            corner[j] = mMin[j] + 0.5 * (mMax[j] - mMin[j]) * (double)(bin_i[j]);
        }

        // The new boxes are half the width in each dimension of this one.
        c_vector<double, DIM> min = corner;
        c_vector<double, DIM> max = corner + 0.5 * (mMax - mMin);

        ParameterBox<DIM>* daughter_box = new ParameterBox<DIM>(this, min, max);

        mDaughterBoxes.push_back(daughter_box);
        CornerSet new_corners_for_this_daughter = daughter_box->GetNewCorners();
        for (unsigned j = 0; j < new_corners_for_this_daughter.size(); j++)
        {
            new_corners.insert(new_corners_for_this_daughter.begin(), new_corners_for_this_daughter.end());
        }
    }

    // Work out interpolated estimates for each QoI
    for (CornerSetIter iter = new_corners.begin();
         iter != new_corners.end();
         ++iter)
    {
        c_vector<double, DIM>* new_corner = *(iter);

        // For each new corner generate an estimate of its QoIs
        // based on interpolation of the existing box.
        std::vector<double> predicted_qois;
        InterpolatePoint(*new_corner, predicted_qois);

        // Store these estimates for comparison with real data later.
        // We don't care whether an error actually occurred or not for this.
        boost::shared_ptr<ParameterPointData> predicted_data = boost::shared_ptr<ParameterPointData>(new ParameterPointData(predicted_qois, 0u));
        AssignQoIValues(new_corner, predicted_data, true);
    }

    // Tidy up things that a parent box doesn't need.
    mCorners.clear();
    mParameterPointDataMapPredictions.clear();
    if (!(mpGreatGrandParentBox == this))
    {
        // Great grandparent box keeps track of the entire data structure,
        // everyone else wipes it on sub-division.
        mParameterPointDataMap.clear();
    }
    mAmParent = true; // I am not actually a parent until I clear my own corners!
    return new_corners;
}

template <unsigned DIM>
void ParameterBox<DIM>::AssignQoIValues(c_vector<double, DIM>* pCorner,
                                        boost::shared_ptr<ParameterPointData> pParameterPointData,
                                        bool isPredictedQoI)
{
    DataMap* p_data_map;
    // If this is a QoI prediction then store it in the predictions map,
    // otherwise store in the 'proper' map.
    if (isPredictedQoI)
    {
        p_data_map = &mParameterPointDataMapPredictions;
    }
    else
    {
        p_data_map = &mParameterPointDataMap;
    }

    // Each corner is only in each box once, so we can safely take the first find result...
    DataMapIter iter = p_data_map->find(pCorner);

    // If we found the relevant corner in this box, and it doesn't already have data associated with it.
    if (iter != p_data_map->end() && (*iter).second == NULL)
    {
        (*iter).second = pParameterPointData;

        // If this is real data and we are a child box, look to see if the fake data exists
        if (!isPredictedQoI && !(mpGreatGrandParentBox == this))
        {
            DataMapIter iter2 = mParameterPointDataMapPredictions.find(pCorner);
            if (iter2 != mParameterPointDataMapPredictions.end())
            {
                // If it does, we found a prediction, so evaluate it.
                std::vector<double> predictions = (*iter2).second->rGetQoIs();
                std::vector<double> real_values = pParameterPointData->rGetQoIs();
                std::vector<double> errors_in_predicitons;
                assert(predictions.size() == real_values.size());
                for (unsigned i = 0; i < predictions.size(); i++)
                {
                    double difference = predictions[i] - real_values[i];
                    //std::cout << "QoI[" << i << "]: prediction = " << predictions[i] << ", \treal = " <<
                    //        real_values[i] << ", \tdifference = " << difference << "\n" << std::flush;
                    errors_in_predicitons.push_back(difference);
                }
                pParameterPointData->SetErrorEstimates(errors_in_predicitons);
                mErrorsInQoIs.push_back(errors_in_predicitons);

                // Erase the predicted ParameterPointData in the map, and remove the entry.
                mParameterPointDataMapPredictions.erase(pCorner);

                // If we have now evaluated all of the predictions
                if (mParameterPointDataMapPredictions.size() == 0)
                {
                    mAllCornersEvaluated = true;

                    mMaxErrorsInEachQoI.clear();
                    // Take each QoI error at first corner to be the max for now.
                    for (unsigned i = 0; i < mErrorsInQoIs[0].size(); i++)
                    {
                        mMaxErrorsInEachQoI.push_back(fabs(mErrorsInQoIs[0][i]));
                    }

                    // For each other corner
                    for (unsigned i = 1u; i < mErrorsInQoIs.size(); i++)
                    {
                        // For each QoI.
                        for (unsigned j = 0; j < mErrorsInQoIs[i].size(); j++)
                        {
                            if (fabs(mErrorsInQoIs[i][j]) > mMaxErrorsInEachQoI[j])
                            {
                                mMaxErrorsInEachQoI[j] = fabs(mErrorsInQoIs[i][j]);
                            }
                        }
                    }
                }
            }
        }
    }

    // Also pass on the instruction to our children boxes.
    for (unsigned i = 0; i < mDaughterBoxes.size(); i++)
    {
        mDaughterBoxes[i]->AssignQoIValues(pCorner, pParameterPointData, isPredictedQoI);
    }
}

template <unsigned DIM>
bool ParameterBox<DIM>::DoesBoxNeedFurtherRefinement(const double& rTolerance,
                                                     const unsigned& rQuantityIndex)
{
    // If I am not the GreatGrandParent
    if (mpParentBox)
        assert(mAllCornersEvaluated);
    if (mpParentBox)
        assert(mMaxErrorsInEachQoI.size() > 0u);
    assert(!mAmParent);
    return (GetMaxErrorInQoIEstimateInThisBox(rQuantityIndex) > rTolerance);
}

template <unsigned DIM>
double ParameterBox<DIM>::GetMaxErrorInQoIEstimateInThisBox(const unsigned& rQuantityIndex)
{
    // If I am the great grandparent box, then we have no idea what the errors are.
    // Could be massive so return DBL_MAX
    if (!mpParentBox)
    {
        return DBL_MAX;
    }

    // Otherwise all these things should be true!
    assert(mAllCornersEvaluated);
    assert(mParameterPointDataMapPredictions.size() == 0u);
    assert(mMaxErrorsInEachQoI.size() > 0u);
    assert(!mAmParent);

    // First check to see if all the corners have errors, if they do we don't want to
    // bother refining this box, so we say it has no error associated with it.
    bool all_errors = true;
    for (DataMapIter iter = mParameterPointDataMap.begin();
         iter != mParameterPointDataMap.end();
         ++iter)
    {
        // See if an error occurred when evaluating the QoIs at this corner.
        if ((*iter).second->GetErrorCode() == 0u) // (Error code 0u means no error occurred)
        {
            all_errors = false;
            break;
        }
    }
    if (all_errors)
    {
        return 0.0;
    }

    // Old code - based on the max gradient of the QoI across the box.
    //    double max = -DBL_MAX;
    //    double min = DBL_MAX;
    //    // This loop will update the maximum and minimum QoI that were recorded, only for corners
    //    // that were not associated with an error in the postprocessing.
    //    for (DataMapIter iter = mParameterPointDataMap.begin();
    //            iter != mParameterPointDataMap.end();
    //            ++iter)
    //    {
    //        const std::vector<double>& r_qois_at_this_parameter_point = (*iter).second->rGetQoIs();
    //        assert(r_qois_at_this_parameter_point.size() >=  quantityIndex + 1u);
    //
    //        if (r_qois_at_this_parameter_point[quantityIndex] > max)
    //        {
    //            max = r_qois_at_this_parameter_point[quantityIndex];
    //        }
    //        if (r_qois_at_this_parameter_point[quantityIndex] < min)
    //        {
    //            min = r_qois_at_this_parameter_point[quantityIndex];
    //        }
    //    }
    //    return max - min;

    // New QoI error-estimate based measure
    return mMaxErrorsInEachQoI[rQuantityIndex];
}

template <unsigned DIM>
void ParameterBox<DIM>::GetErrorEstimateInAllBoxes(ParameterBox<DIM>*& pBestBox,
                                                   double& rErrorEstimateInBestBox,
                                                   const double& rTolerance,
                                                   const unsigned& rQuantityIndex)
{
    if (!mAmParent)
    {
        // If I have a parent (if I am not the great grandparent)
        if (mpParentBox)
            assert(mAllCornersEvaluated);
        if (mpParentBox)
            assert(mMaxErrorsInEachQoI.size() > 0u);
        assert(mParameterPointDataMap.size() == pow(2, DIM)); // Check we have an entry for all of our corners.
        assert(mDaughterBoxes.size() == 0u);

        if (DoesBoxNeedFurtherRefinement(rTolerance, rQuantityIndex))
        {
            double max_error_estimate = GetMaxErrorInQoIEstimateInThisBox(rQuantityIndex);
            unsigned errors_associated = GetNumErrors();

            // If we have a new maximum variation and
            //   * we have no current best box,
            //   * (there may be an error but we'll have to go with it.)
            // OR we have a new maximum variation and
            //   * a current best box guess,
            //   * we have more than one error in the evaluation of the results in this box.
            bool new_max = (max_error_estimate > rErrorEstimateInBestBox);
            if ((new_max && !pBestBox) || // This condition applies to the great-grandparent original box.
                (new_max && errors_associated <= pow(2, DIM) - 1)) // || // We don't want to start refining where there are loads of errors
            //(errors_associated < pBestBox->GetNumErrors()) ) // So boxes with no error always win
            {
                pBestBox = this;
                rErrorEstimateInBestBox = max_error_estimate;
            }
        }
    }
    else
    {
        // If I am a parent then ask my daughter boxes the same thing.
        for (unsigned i = 0; i < mDaughterBoxes.size(); i++)
        {
            mDaughterBoxes[i]->GetErrorEstimateInAllBoxes(pBestBox, rErrorEstimateInBestBox, rTolerance, rQuantityIndex);
        }
    }
    return;
}

template <unsigned DIM>
ParameterBox<DIM>* ParameterBox<DIM>::FindBoxWithLargestQoIErrorEstimate(const unsigned& rQuantityIndex,
                                                                         const double& rTolerance,
                                                                         const unsigned& rMaxGenerationDifference)
{
    // Only the grand parent should call this. If I have a parent I'm not it.
    if (mpParentBox)
    {
        EXCEPTION("Only the original parameter box should call this method.");
    }

    ParameterBox<DIM>* p_box = NULL;
    double variation = -DBL_MAX;
    GetErrorEstimateInAllBoxes(p_box, variation, rTolerance, rQuantityIndex);

    // If there is somewhere that doesn't meet the tolerances
    if (variation > rTolerance)
    {
        assert(p_box);

        ParameterBox<DIM>* most_refined = GetMostRefinedChild();
        ParameterBox<DIM>* least_refined = GetLeastRefinedChild(rTolerance, rQuantityIndex);

        // Check the selected box isn't going to refine one area too much,
        // if it is refine least refined area instead.
        if (least_refined // if an unrefined box exists that doesn't meet the tolerances.
            && ((most_refined->GetGeneration() - least_refined->GetGeneration()) == rMaxGenerationDifference)
            && (p_box->GetGeneration() == GetMostRefinedChild()->GetGeneration()))
        {
            return least_refined;
        }
    }

    // otherwise
    return p_box;
}

template <unsigned DIM>
ParameterBox<DIM>* ParameterBox<DIM>::GetBoxContainingPoint(const c_vector<double, DIM>& rPoint)
{
    // We only want to perform this check on the main box.
    if (!mpParentBox && !this->IsPointInThisBox(rPoint))
    {
        EXCEPTION("This point is not contained within this box (or any of its children).");
    }

    // This box definitely contains the point
    // (if we aren't the original box, this was already checked by this method lower down before recursive call!)
    // If it isn't a parent then it is the box we are looking for.
    if (!mAmParent)
    {
        assert(this->IsPointInThisBox(rPoint));
        return this;
    }

    // If it is a parent then ask daughters...
    ParameterBox<DIM>* p_box = this;

    std::vector<ParameterBox<DIM>*> daughters = GetDaughterBoxes();
    for (unsigned i = 0; i < daughters.size(); i++)
    {
        if (daughters[i]->IsPointInThisBox(rPoint))
        {
            // Recursively call this method on each daughter box.
            p_box = daughters[i]->GetBoxContainingPoint(rPoint);
        }
    }
    return p_box;
}

template <unsigned DIM>
std::vector<double> ParameterBox<DIM>::InterpolateQoIsAt(const c_vector<double, DIM>& rPoint)
{
    // Only the grand parent should call this.
    if (mpParentBox)
    {
        EXCEPTION("Only the original parameter box should call this method.");
    }

    std::vector<double> interpolated_qois;

    ParameterBox<DIM>* p_box = GetBoxContainingPoint(rPoint);
    p_box->InterpolatePoint(rPoint, interpolated_qois);

    return interpolated_qois;
}

template <unsigned DIM>
void ParameterBox<DIM>::InterpolatePoint(const c_vector<double, DIM>& rPoint,
                                         std::vector<double>& rQoIs)
{
    // I am a child box.
    assert(!mAmParent);
    // And check some data has been assigned to my corners
    assert((*(mParameterPointDataMap.begin())).second != boost::shared_ptr<ParameterPointData>());

    // Work out how many QoIs we should have, same across all parameter points, so just take the first one.
    DataMapIter iter = mParameterPointDataMap.begin();
    unsigned num_qois = (*iter).second->rGetQoIs().size();

    // Wipe the existing QoIs vector.
    rQoIs.resize(num_qois);
    for (unsigned qoi_idx = 0; qoi_idx < num_qois; qoi_idx++)
    {
        rQoIs[qoi_idx] = 0.0;
    }

    c_vector<double, DIM> point = rPoint;
    // Nondimensionalise the point within this box
    for (unsigned j = 0; j < DIM; j++)
    {
        point[j] = (rPoint[j] - mMin[j]) / (mMax[j] - mMin[j]);
    }

    // See doxygen comment for this method for some detail of what is going on here.
    // This is in the same order as corners (as it is how we originally calculated them!).
    for (unsigned i = 0; i < pow(2, DIM); i++)
    {
        // Use a binary conversion to get the right indices in place.
        std::bitset<DIM> bin_i(i);
        double multiplier_for_this_corner = 1.0;
        for (unsigned j = 0; j < DIM; j++)
        {
            if (bin_i[j] == 0)
            {
                multiplier_for_this_corner *= (1.0 - point[j]);
            }
            else
            {
                multiplier_for_this_corner *= point[j];
            }
        }

        for (unsigned qoi_idx = 0; qoi_idx < num_qois; qoi_idx++)
        {
            rQoIs[qoi_idx] += multiplier_for_this_corner * mParameterPointDataMap[mCorners[i]]->rGetQoIs()[qoi_idx];
        }
    }
}

template <unsigned DIM>
bool ParameterBox<DIM>::IsPointInThisBox(const c_vector<double, DIM>& rPoint)
{
    bool within = true;
    for (unsigned i = 0; i < DIM; i++)
    {
        if (rPoint[i] > mMax[i] || rPoint[i] < mMin[i])
        {
            within = false;
            break;
        }
    }
    return within;
}

template <unsigned DIM>
unsigned ParameterBox<DIM>::GetNumErrors()
{
    unsigned error_count = 0u;

    for (unsigned i = 0; i < mCorners.size(); i++)
    {
        if (mParameterPointDataMap[mCorners[i]]->GetErrorCode() > 0u)
        {
            error_count++;
        }
    }

    return error_count;
}

template <unsigned DIM>
ParameterBox<DIM>* ParameterBox<DIM>::GetMostRefinedChild()
{
    ParameterBox<DIM>* p_box = NULL;

    if (!mAmParent)
    {
        assert(mCorners.size() == pow(2, DIM));
        assert(mDaughterBoxes.size() == 0u);
        p_box = this;
    }
    else
    {
        // Get a load of the most refined boxes that are owned by daughter boxes.
        for (unsigned i = 0; i < mDaughterBoxes.size(); i++)
        {
            ParameterBox<DIM>* this_daughters_most_refined = mDaughterBoxes[i]->GetMostRefinedChild();
            if (!p_box || this_daughters_most_refined->GetGeneration() > p_box->GetGeneration())
            {
                p_box = this_daughters_most_refined;
            }
        }
    }

    return p_box;
}

template <unsigned DIM>
ParameterBox<DIM>* ParameterBox<DIM>::GetLeastRefinedChild(const double& rTolerance, const unsigned& rQuantityIndex)
{
    ParameterBox<DIM>* p_box = NULL;

    if (!mAmParent)
    {
        assert(mCorners.size() == pow(2, DIM));
        assert(mDaughterBoxes.size() == 0u);

        if (DoesBoxNeedFurtherRefinement(rTolerance, rQuantityIndex))
        {
            p_box = this;
        }
    }
    else
    {
        // Get a load of the least refined boxes that are owned by daughter boxes.
        for (unsigned i = 0; i < mDaughterBoxes.size(); i++)
        {
            ParameterBox<DIM>* this_daughters_least_refined = mDaughterBoxes[i]->GetLeastRefinedChild(rTolerance, rQuantityIndex);

            if (this_daughters_least_refined // If this daughter box needs refinement
                && (!p_box // and we either don't have a box at the moment, or this box is a lower generation
                    || this_daughters_least_refined->GetGeneration() < p_box->GetGeneration()))
            {
                p_box = this_daughters_least_refined;
            }
        }
    }

    return p_box;
}

template <unsigned DIM>
std::vector<double> ParameterBox<DIM>::GetMaxErrorsInPredictedQoIs() const
{
    if (!mAllCornersEvaluated)
    {
        EXCEPTION("Not all the parameter points (which you can get with GetNewCorners()) have been assigned data. Error estimates unavailable.");
    }

    assert(mMaxErrorsInEachQoI.size() > 0u);
    assert(mParameterPointDataMapPredictions.size() == 0u);
    return mMaxErrorsInEachQoI;
}

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ParameterBox)

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class ParameterBox<1>;
template class ParameterBox<2>;
template class ParameterBox<3>;
template class ParameterBox<4>;
template class ParameterBox<5>;
