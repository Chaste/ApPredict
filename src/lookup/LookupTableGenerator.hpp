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

#ifndef LOOKUPTABLEGENERATOR_HPP_
#define LOOKUPTABLEGENERATOR_HPP_

#include <set>
#include <boost/shared_ptr.hpp>
// Seems that whatever version of ublas I am using now contains boost serialization
// methods for c_vector, which is nice.
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "ChasteSerialization.hpp" // Should be included before any other Chaste headers.
#include "ChasteSerializationVersion.hpp"

#include "AbstractCvodeCell.hpp"
#include "ParameterBox.hpp"
#include "Exception.hpp"
#include "UblasVectorInclude.hpp" // Chaste helper header to get c_vectors included with right namespace.
#include "OutputFileHandler.hpp"
#include "SingleActionPotentialPrediction.hpp"
#include "QuantityOfInterest.hpp"
#include "ParameterPointData.hpp"

/**
 * A class that will generate lookup tables in DIM-dimensional parameter space,
 * populate them, perform refinement, etc.
 *
 * Originally written to make a lookup table for percentage block of
 *
 * IKr  (hERG)
 * INa  (NaV1.5)
 * ICaL (CaV1.2)
 * IKs  (KCNQ1/minK)
 * Ito  (Kv4.3/KChIP2.2)
 *
 * You should add QoIs in order of importance, as the lookup table will be refined for each in turn.
 */
template<unsigned DIM>
class LookupTableGenerator
{
  private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mModelIndex;
        if (version > 0u)
        {
            archive & mFrequency;
        }
        archive & mParameterPoints;
        archive & mParameterPointData;
        archive & mParameterNames;
        archive & mMinimums;
        archive & mMaximums;
        archive & mUnscaledParameters;
        archive & mQuantitiesToRecord;
        archive & mQoITolerances;
        archive & mMaxNumEvaluations;
        archive & mNumEvaluations;
        archive & mOutputFileName;
        archive & mOutputFolder;
        archive & mGenerationHasBegun;
        archive & mpParentBox;
        archive & mInitialConditions;
        archive & mMaxRefinementDifference;
        archive & mMaxNumPaces;
    }

    /** Helper wrappers round these long-winded set and iterator names */
    typedef typename std::set<c_vector<double, DIM>*, c_vector_compare<DIM> > CornerSet;
    typedef typename std::set<c_vector<double, DIM>*, c_vector_compare<DIM> >::iterator CornerSetIter;

    /** The model index we are using, refers to options in SetupModel class. */
    unsigned mModelIndex;

    /** The frequency of pacing to use (in Hz), default value is taken to be 1Hz */
    double mFrequency;

    /** The steady 1Hz pacing conditions for ones(DIM). */
    std::vector<double> mInitialConditions;

    /** The points in DIM-dimensional parameter space at which we have evaluated quantities of interest (lookup table co-ords) */
    std::vector<c_vector<double, DIM> > mParameterPoints;

    /** The data (QoIs) that were evaluated at the points in #mParameterPoints */
    std::vector<boost::shared_ptr<ParameterPointData> > mParameterPointData;

    /** The Oxford metadata names of the parameters that we are going to scale */
    std::vector<std::string> mParameterNames;

    /** The minimum scalings that should be applied (lower bounds of lookup table) */
    std::vector<double> mMinimums;

    /** The maximum scalings that should be applied (upper bounds of lookup table) */
    std::vector<double> mMaximums;

    /** The unscaled parameters that should not be messed with! */
    std::vector<double> mUnscaledParameters;

    /** The list of quantities of interest that we should record */
    std::vector<QuantityOfInterest> mQuantitiesToRecord;

    /** The tolerances to aim for for each Quantity of Interest (QoI). */
    std::vector<double> mQoITolerances;

    /** The maximum number of function evaluations to perform */
    unsigned mMaxNumEvaluations;

    /** The number of evaluations that have been performed up to now. */
    unsigned mNumEvaluations;

    /** The output file base name to use */
    std::string mOutputFileName;

    /** The output folder to use when writing out lookup table entries. */
    std::string mOutputFolder;

    /** Whether we have started generating lookup tables, to prevent accidents */
    bool mGenerationHasBegun;

    /** The maximum difference in parameter box refinement we are going to allow */
    unsigned mMaxRefinementDifference;

    /** The great-great-grandaddy parameter box */
    ParameterBox<DIM>* mpParentBox;

    /** The maximum number of paces to do */
    unsigned mMaxNumPaces;

    /**
     * This method will farm out the evaluation of a set of points using multi-threading.
     *
     * @param setOfPoints  A collection of points in parameter space at which to evaluate QoIs.
     * @param rFile  The output file to write a line of results into.
     */
    void RunEvaluationsForThesePoints(CornerSet setOfPoints, out_stream& rFile);

    /**
     *  Private constructor, just for use in archiving
     *
     *  This just ensures the pointers are correctly initialised to
     *  make sure we can detect and recover the intended state of everything.
     */
    LookupTableGenerator()
     : mModelIndex(0u),
       mpParentBox(NULL)
    {};

  public:

    /**
     * Constructor
     *
     * @param rModelIndex  The action potential model to use (see SetupModel for options).
     * @param rOutputFileName  The base file name to use for the lookup table.
     * @param rOutputFolder  The folder to put the output in, defaults to "LookupTables".
     */
    LookupTableGenerator(const unsigned& rModelIndex,
                         const std::string& rOutputFileName,
                         const std::string& rOutputFolder = "LookupTables");

    /**
     * Destructor - cleans up some memory
     */
    ~LookupTableGenerator();

    /**
     * Main method,
     *
     * runs simulations varying the parameters requested in SetParameterToScale()
     *
     * performs postprocessing requested by AddQuantityOfInterest()
     *
     * and puts lookup table entries for f(mParameterPoints) =  mQuantitiesOfInterest
     * in the relevant member variables.
     */
    void GenerateLookupTable();

    /**
     * @return The points in NUM_PARAMS-dimensional parameter space at which the quantities of interest have been evaluated.
     */
    std::vector<c_vector<double, DIM> > GetParameterPoints();

    /**
     * @return The quantities of interest that have been evaluated at points in the lookup table.
     */
    std::vector<std::vector<double> > GetFunctionValues();

    /**
     * Choose which parameter in the model to scale.
     *
     * @param rMetadataName  The Oxford metadata name of the parameter we want to scale.
     * @param rMin  The minimum bound the parameter should take (limit of lookup table).
     * @param rMax  The maximum bound the parameter should take (limit of lookup table).
     */
    void SetParameterToScale(const std::string& rMetadataName,
                             const double& rMin,
                             const double& rMax);

    /**
     * Set the maximum number of paces to do in a single evaluation.
     *
     * This changes the "timeout" on the steady state calculation. If we
     * reach steady state first then we stop there.
     *
     * This method is useful if you are emulating an experiment of a given
     * duration. For example in the rabbit wedge at GSK they don't go for
     * longer than 30 minutes at a given frequency so we change this accordingly.
     *
     * @param numPaces  The maximum number of paces to do.
     */
    void SetMaxNumPaces(unsigned numPaces);

    /**
     * @return The maximum number of paces considered in the Lookup table generation.
     */
    unsigned GetMaxNumPaces();

    /**
     * Add a quantity of interest to create a lookup table for.
     *
     * These are listed in QuantityOfInterest enumeration at the top of LookupTableGenerator.hpp
     *
     * @param quantity  A quantity of interest (QoI)
     * @param tolerance the tolerance we are willing to accept in this QoI between lookup table points.
     */
    void AddQuantityOfInterest(QuantityOfInterest quantity, double tolerance);

    /**
     * Set maximum number of simulations to run when making the table.
     *
     * @param rMaxNumEvals  The maximum number of simulations.
     */
    void SetMaxNumEvaluations(const unsigned& rMaxNumEvals);

    /**
     * Set a maximum difference in refinement level across the parameter space.
     *
     * If you don't set this there is assumed to be no limit on the refinement around the areas
     * where the QoIs are varying most rapidly, this can lead to under-exploring other parts of
     * the parameter space.
     *
     * @param rMaxRefinementDifference  The maximum difference in the parameter boxes to allow.
     */
    void SetMaxVariationInRefinement(const unsigned& rMaxRefinementDifference);

    /**
     * Provide an interpolated estimate for the quantities of interest throughout parameter space.
     *
     * @param rParameterPoints  The points in parameter space at which we would like to estimate QoIs.
     * @return The QoI estimates at these points.
     */
    std::vector<std::vector<double> > Interpolate(const std::vector<c_vector<double, DIM> >& rParameterPoints);

    /**
     * @return The number of evaluations (points in the lookup table at which Quantities of Interest have been evaluated).
     */
    unsigned GetNumEvaluations();

    /**
     * Set the pacing frequency to use throughout the lookup table generation.
     *
     * @param frequency  The pacing frequency to use (in Hz).
     */
    void SetPacingFrequency(double frequency);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LookupTableGenerator)

// Keep track of the archive version we are using
namespace boost
{
    namespace serialization
    {
        /**
         * Specify a version number for archive backwards compatibility.
         *
         * This is how to do BOOST_CLASS_VERSION(AbstractCardiacPde, 1)
         * with a templated class.
         */
        template <unsigned DIM>
        struct version<LookupTableGenerator<DIM> >
        {
            CHASTE_VERSION_CONTENT(1); // Increment this on serialize method changes.
        };
    } // namespace serialization
} // namespace boost


#endif // LOOKUPTABLEGENERATOR_HPP_
