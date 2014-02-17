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

#ifndef PARAMETERBOX_HPP
#define PARAMETERBOX_HPP

#include <set>
#include <map>
#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"  // Should be included before any other Chaste headers.
#include "ParameterPointData.hpp"
#include "UblasVectorInclude.hpp" // Chaste helper header to get c_vectors included with right namespace.

// Seems that whatever version of ublas I am using now contains boost serialization
// methods for c_vector, which is nice.
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

static const double TOL = 1e-12;

/**
 *  A special comparison method to allow std::map to sort and compare
 *  a std::map<c_vector<double,DIM>* >
 */
template<unsigned DIM>
struct c_vector_compare
{
    // Sorts on x, then y, then ... etc.
    bool operator () (const c_vector<double,DIM>* a, const c_vector<double,DIM>* b) const
    {
        bool a_is_smaller = false;
        for (unsigned i=0; i<DIM; i++)
        {
            if ((*a)[i] < (*b)[i] - TOL)
            {
                a_is_smaller = true;
                break;
            }
            if ((*a)[i] > (*b)[i] + TOL)
            {
                a_is_smaller = false;
                break;
            }
            // We only carry on to next index if they are equal.
        }
        return a_is_smaller;
    }
};

/**
 * This class stores the co-ordinates of the corners of N-D boxes.
 *
 * It provides a method for dividing a box into equally sized sub-boxes, *
 * and keeps track of their indices too.
 *
 * It also provides methods to find out which sub-boxes have the most and
 * least refinement, as well as suggesting the next one to refine.
 */
template<unsigned DIM>
class ParameterBox
{
  private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    friend class TestParameterBox;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mAmParent;
        archive & mpParentBox;
        archive & mpGreatGrandParentBox;
        archive & mMin;
        archive & mMax;
        archive & mCorners;
        archive & mNewCorners;
        archive & mDaughterBoxes;
        archive & mParameterPointDataMap;
        archive & mParameterPointDataMapPredictions;
        archive & mErrorsInQoIs;
        archive & mMaxErrorsInEachQoI;
        archive & mAllCornersEvaluated;
        archive & mGeneration;
    }

    typedef typename std::set<c_vector<double, DIM>*, c_vector_compare<DIM> > CornerSet;
    typedef typename std::set<c_vector<double, DIM>*, c_vector_compare<DIM> >::iterator CornerSetIter;
    typedef typename std::vector<c_vector<double, DIM>* > CornerVec;
    typedef typename std::vector<c_vector<double, DIM>* >::iterator CornerVecIter;
    typedef typename std::map<c_vector<double, DIM>*, boost::shared_ptr<ParameterPointData>, c_vector_compare<DIM> > DataMap;
    typedef typename std::map<c_vector<double, DIM>*, boost::shared_ptr<ParameterPointData>, c_vector_compare<DIM> >::iterator DataMapIter;

    /** Whether this box has been subdivided and has children boxes */
    bool mAmParent;

    /** This box's parent box */
    ParameterBox<DIM>* mpParentBox;

    /** The original parent box */
    ParameterBox<DIM>* mpGreatGrandParentBox;

    /** The number of 'subdivisions' we have done to get this box. Top box is generation 0.*/
    unsigned mGeneration;

    /** The N-Dim minimum of this box */
    c_vector<double, DIM> mMin;
    /** The N-Dim maximum of this box */
    c_vector<double, DIM> mMax;

    /** The locations of this box's vertices in N-D space */
    std::vector<c_vector<double, DIM>* > mCorners;

    /** New corners that were created especially for this box */
    CornerSet mNewCorners;

    /** Pointers to the children (subdivisions) of this box */
    std::vector<ParameterBox<DIM>* > mDaughterBoxes;

    /**
     * The QoIs and errors associated with this box's corners,
     * only populated by the great-grandparent and the children
     * not the intermediate boxes.
     *
     * Children keep a minimal map of their own corners,
     *
     * Great-grandparent keeps an entire map of the whole thing.
     */
    DataMap mParameterPointDataMap;

    /**
     * A map of predicted ParameterPointData for new
     * corners, populated at subdivision time.
     *
     * This is only held by the local children boxes.
     */
    DataMap mParameterPointDataMapPredictions;

    /**
     * The errors in the QoIs at each new corner / predicted data point.
     * This is unordered (but the values are propagated into the relevant
     * ParameterPointData structure).
     *
     * It is populated each time you call AssignQoI values.
     */
    std::vector<std::vector<double> > mErrorsInQoIs;

    /**
     * The maximum errors in each QoIs over all new data points.
     *
     * This is populated by AssignQoI values when it has been called for
     * ALL of the new corners.
     */
    std::vector<double> mMaxErrorsInEachQoI;

    /**
     * Whether the error estimates are complete (all corners evaluated).
     */
    bool mAllCornersEvaluated;

    /**
     * Get a measure of the error associated with predicting the quantity of interest in this box.
     * Calculated by providing an interpolated estimate from parent, and then comparing with
     * our new assigned quantities of interest.
     *
     * Any QoI that were associated with errors on evaluation are ignored,
     * and error is said to be zero (no refinement needed).
     *
     * @param quantityIndex  The quantity of interest.
     * @return The largest difference between interpolated and evaluated that this quantity takes in this box.
     */
    double GetMaxErrorInQoIEstimateInThisBox(const unsigned& rQuantityIndex);

    /**
     * Calls GetMaxErrorInQoIEstimateInThisBox() for any children and updates a pointer to the
     * one with the largest error estimate above the tolerance.
     *
     * @param pBestBox  Reference to a pointer to the box with the largest variation at present.
     * @param rErrorEstimateInBestBox  The value of the largest error estimate associated with the best box.
     * @param rTolerance  The error estimate we are happy with.
     * @param rQuantityIndex  The index of the quantity of interest we are examining at present.
     */
    void GetErrorEstimateInAllBoxes(ParameterBox<DIM>*& pBestBox,
                                    double& rErrorEstimateInBestBox,
                                    const double& rTolerance,
                                    const unsigned& rQuantityIndex);

    /** Private constructor, just for use by archiving */
    ParameterBox(){};

    /**
     * @return A pointer to (one of) the children with the largest generation number.
     */
    ParameterBox<DIM>* GetMostRefinedChild();

    /**
     * Get the least refined child box that doesn't meet the tolerances. If all meet the tolerances it
     * will return NULL.
     *
     * @param rTolerance  The error estimate we are happy with.
     * @param rQuantityIndex  The index of the quantity of interest we are examining at present.
     * @return A pointer to (one of) the children (that has no children of its own)
     *         with the smallest generation number.
     */
    ParameterBox<DIM>* GetLeastRefinedChild(const double& rTolerance, const unsigned& rQuantityIndex);

    /**
     * @return Whether a point is contained within a box, or one of its children.
     */
    bool IsPointInThisBox(const c_vector<double, DIM>& rPoint);

    /**
     * Perform interpolation within this box.
     *
     * The algorithm implements something like the following example for 3D:
     * taken from http://paulbourke.net/miscellaneous/interpolation/
     *
     * Given a non-dimensionalised point x,y,z, in an N-DIM box of sides 0->1.
     *
     * V_xyz =   V_000 (1 - x) (1 - y) (1 - z) +
     *           V_100 x       (1 - y) (1 - z) +
     *           V_010 (1 - x) y       (1 - z) +
     *           V_001 (1 - x) (1 - y) z       +
     *           V_101 x       (1 - y) z       +
     *           V_011 (1 - x) y       z       +
     *           V_110 x       y       (1 - z) +
     *           V_111 x       y       z
     *
     * This fits in really nicely with the way we specified the corners in the first place
     * with binary operations (see constructor). This has allowed us to generalise to N-DIM
     * really easily.
     *
     * @param rPoint  a point in parameter space at which to interpolate within THIS box.
     * @param rQoIs  a reference to QoIs, populated by this method.
     */
    void InterpolatePoint(const c_vector<double, DIM>& rPoint, std::vector<double>& rQoIs);

    /**
     * @return The generation of this box.
     */
    unsigned GetGeneration();

    /**
     * @return The number of errors associated with evaluating the corners of this box.
     */
    unsigned GetNumErrors();

    /**
     * Get all the parameter points contained within this box and all its children as a vector
     *
     * @return points in DIM-dimensional parameter space
     */
    std::vector<c_vector<double, DIM>* > GetCornersAsVector();

    /**
     * Get all the parameter points contained at corners of this box
     *
     * @return points in DIM-dimensional parameter space
     */
    std::vector<c_vector<double, DIM>* > GetOwnCorners();

    /**
     * Get all the parameter points created especially for this box
     *
     * @return points in DIM-dimensional parameter space
     */
    CornerSet GetNewCorners();

    /**
     * @return Whether this box has been subdivided and has children boxes.
     */
    bool IsParent();

    /**
     * @return This box's parent box
     */
    ParameterBox<DIM>* GetParent();

    /**
     * @return This box's original great-grand parent box
     */
    ParameterBox<DIM>* GetGreatGrandParent();

    /**
     * @return The daughter boxes that are contained (directly) within this box.
     *
     * Does not include any grand children etc.
     */
    std::vector<ParameterBox<DIM>* > GetDaughterBoxes();

    /**
     * @return All the boxes that are contained within this box.
     *
     * Includes all grand children etc.
     */
    std::vector<ParameterBox<DIM>* > GetWholeFamilyOfBoxes();

    /**
     * Find the parameter box that contains a certain point in DIM-dimensional parameter space.
     *
     * If no boxes containing the point exist then an exception is thrown.
     *
     * @param rPoint  The point of interest.
     *
     * @return The box that contains the point of interest.
     */
    ParameterBox<DIM>* GetBoxContainingPoint(const c_vector<double, DIM>& rPoint);

    /**
     * Whether this box needs further refinement or not.
     *
     * @param rTolerance  The error estimate in the QoI that we are happy with
     * @param rQuantityIndex  The QoI index that we are considering.
     * @return Whether this box needs further refinement or not.
     */
    bool DoesBoxNeedFurtherRefinement(const double& rTolerance, const unsigned& rQuantityIndex);

  public:
    /**
     * Constructor
     *
     * @param pParent The parent of this box
     * @param rMin  The N-D minimum vertex of this box
     * @param rMax  THe N-D maximum vertex of this box
     */
    ParameterBox(ParameterBox<DIM>* pParent,
                 const c_vector<double, DIM>& rMin = zero_vector<double>(DIM),
                 const c_vector<double, DIM>& rMax = scalar_vector<double>(DIM, 1.0));

    /**
     * Destructor - deletes children boxes
     */
    ~ParameterBox();

    /**
     * Get all the parameter points contained within this box and all its children.
     *
     * @return points in DIM-dimensional parameter space
     */
    CornerSet GetCorners();

    /**
     * Subdivide this box into 2^DIM, e.g.
     *  * 1D - two
     *  * 2D - four
     *  * 3D - eight
     *
     * Makes this box into a parent and creates daughter boxes.
     *
     * @return The new points in parameter space at which quantities of interest need to be evaluated.
     */
    CornerSet SubDivide();


    /**
     * Tell this box, and any children, the values of the quantities of interest (QoIs) at a given point.
     *
     * @param pCorner  The point in parameter space at which we have evaluated the QoIs.
     * @param pParameterPointData The data values at this point in parameter space.
     * @param isPredictedQoI  Whether this is a predicted QoI or a real one, data is entered into the predicted or real map accordingly.
     */
    void AssignQoIValues(c_vector<double, DIM>* pCorner,
                         boost::shared_ptr<ParameterPointData> pParameterPointData,
                         bool isPredictedQoI = false);

    /**
     * Find the parameter box that has the largest error estimate in a given quantity of interest.
     * Will not over-refine one area if a rMaxGenerationDifference is set.
     *
     * If all boxes meet the tolerance a null pointer is returned.
     *
     * @param quantityIndex  The index of the quantity of interest to check.
     * @param tolerance  The tolerance for this quantity of interest across the box.
     * @param maxGenerationDifference  The maximum difference in refinement levels in terms of generation.
     *
     * @return The box that most exceeds the tolerance in the quantity of interest, if any.
     */
    ParameterBox<DIM>* FindBoxWithLargestQoIErrorEstimate(const unsigned& rQuantityIndex,
                                                          const double& rTolerance,
                                                          const unsigned& rMaxGenerationDifference = UNSIGNED_UNSET);

    /**
     * Calculate a regular grid interpolation by finding the box containing the point and
     * interpolating QoIs from its corners.
     *
     * @param rPoint  The point in DIM-dimensional parameter space at which we want an estimate of the QoI
     * @return  The estimate of the QoIs at this point.
     */
    std::vector<double> InterpolateQoIsAt(const c_vector<double, DIM>& rPoint);

    /**
     * Return the size of the maximum (across all corners) error in predictions of each QoI.
     *
     * @return maximum error in each QoI.
     */
    std::vector<double> GetMaxErrorsInPredictedQoIs() const;
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ParameterBox)



#endif // PARAMETERBOX_HPP
