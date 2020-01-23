/*

Copyright (c) 2005-2020, University of Oxford.
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

#ifndef ABSTRACTUNTEMPLATEDLOOKUPTABLEGENERATOR_HPP_
#define ABSTRACTUNTEMPLATEDLOOKUPTABLEGENERATOR_HPP_

#include <boost/shared_ptr.hpp>

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "AbstractCvodeCell.hpp"
#include "QuantityOfInterest.hpp"
#include "SingleActionPotentialPrediction.hpp"

/**
 * A class to allow a nice interface for ApPredictMethods to use any dimension
 * Lookup table.
 *
 * No member variables...
 */
class AbstractUntemplatedLookupTableGenerator
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Main boost serialization method, some member variables are handled by
     * Save/Load constructs.
     *
     * @param archive  the archive file.
     * @param version  the version of archiving, defined in the macro at the
     * bottom of the class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        // No private members to archive,
        // but helpful to have in hierarchy as we want to archive via a pointer to
        // this class.
    }

public:
    /**
    * Constructor - default.
    */
    AbstractUntemplatedLookupTableGenerator(){};

    /**
     * Destructor - default
     */
    virtual ~AbstractUntemplatedLookupTableGenerator(){};

    /**
      * Main method,
      *
      * runs simulations varying the parameters requested in
      * SetParameterToScale()
      *
      * performs postprocessing requested by AddQuantityOfInterest()
      *
      * and puts lookup table entries for f(mParameterPoints) =
      * mQuantitiesOfInterest
      * in the relevant member variables.
      *
      * @return whether the table meets all tolerances.
      */
    virtual bool GenerateLookupTable() = 0;

    //    /**
    //	 * @return The points in NUM_PARAMS-dimensional parameter space at which
    // the
    //	 * quantities of interest have been evaluated.
    //	 */
    //    template <unsigned DIM>
    //    std::vector<c_vector<double, DIM> > GetParameterPoints();

    /**
      * @return The quantities of interest that have been evaluated at points
      * in
      * the lookup table.
      */
    virtual std::vector<std::vector<double> > GetFunctionValues() = 0;

    /**
      * Choose which parameter in the model to scale.
      *
      * @param rMetadataName  The Oxford metadata name of the parameter we want
      * to
      * scale.
      * @param rMin  The minimum bound the parameter should take (limit of
      * lookup
      * table).
      * @param rMax  The maximum bound the parameter should take (limit of
      * lookup
      * table).
      */
    virtual void SetParameterToScale(const std::string& rMetadataName,
                                     const double& rMin, const double& rMax)
        = 0;

    /**
      * Set the maximum number of paces to do in a single evaluation.
      *
      * This changes the "timeout" on the steady state calculation. If we
      * reach steady state first then we stop there.
      *
      * This method is useful if you are emulating an experiment of a given
      * duration. For example in the rabbit wedge at GSK they don't go for
      * longer than 30 minutes at a given frequency so we change this
      * accordingly.
      *
      * @param numPaces  The maximum number of paces to do.
      */
    virtual void SetMaxNumPaces(unsigned numPaces) = 0;

    /**
      * @return The maximum number of paces considered in the Lookup table
      * generation.
      */
    virtual unsigned GetMaxNumPaces() = 0;

    /**
      * Add a quantity of interest to create a lookup table for.
      *
      * These are listed in QuantityOfInterest enumeration at the top of
      * LookupTableGenerator.hpp
      *
      * @param quantity  A quantity of interest (QoI)
      * @param tolerance the tolerance we are willing to accept in this QoI
      * between
      * lookup table points.
      */
    virtual void AddQuantityOfInterest(QuantityOfInterest quantity, double tolerance) = 0;

    /**
      * Set maximum number of simulations to run when making the table.
      *
      * @param rMaxNumEvals  The maximum number of simulations.
      */
    virtual void SetMaxNumEvaluations(const unsigned& rMaxNumEvals) = 0;

    /**
      * Set a maximum difference in refinement level across the parameter
      * space.
      *
      * If you don't set this there is assumed to be no limit on the
      * refinement
      * around the areas
      * where the QoIs are varying most rapidly, this can lead to
      * under-exploring
      * other parts of
      * the parameter space.
      *
      * @param rMaxRefinementDifference  The maximum difference in the
      * parameter
      * boxes to allow.
      */
    virtual void SetMaxVariationInRefinement(const unsigned& rMaxRefinementDifference) = 0;

    //    /**
    //	 * Provide an interpolated estimate for the quantities of interest
    // throughout
    //	 * parameter space.
    //	 *
    //	 * @param rParameterPoints  The points in parameter space at which we
    // would
    //	 * like to estimate QoIs.
    //	 * @return The QoI estimates at these points.
    //	 */
    //    template <unsigned DIM>
    //    std::vector<std::vector<double> > Interpolate(const
    //    std::vector<c_vector<double, DIM> >& rParameterPoints);

    /**
     * Provide an interpolated estimate for the quantities of interest
     * throughout
     * parameter space.
     *
     * @param rParameterPoints  The points in parameter space at which we
     * would like to estimate QoIs.
     * @return The QoI estimates at these points.
     */
    virtual std::vector<std::vector<double> > Interpolate(const std::vector<std::vector<double> >& rParameterPoints) = 0;

    /**
     * @return The number of evaluations (points in the lookup table at which
     * Quantities of Interest have been evaluated).
     */
    virtual unsigned GetNumEvaluations() = 0;

    /**
     * Set the pacing frequency to use throughout the lookup table generation.
     *
     * @param frequency  The pacing frequency to use (in Hz).
     */
    virtual void SetPacingFrequency(double frequency) = 0;

    /**
     * Helper method that just returns DIM, to avoid template chaos.
     */
    virtual unsigned GetDimension() const = 0;

    /**
     * @return The names of each parameter in each dimension of this table.
     */
    virtual std::vector<std::string> GetParameterNames() const = 0;
};

CLASS_IS_ABSTRACT(AbstractUntemplatedLookupTableGenerator)

#endif // ABSTRACTUNTEMPLATEDLOOKUPTABLEGENERATOR_HPP_
