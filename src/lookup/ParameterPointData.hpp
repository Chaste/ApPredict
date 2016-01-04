/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef PARAMETERPOINTDATA_HPP_
#define PARAMETERPOINTDATA_HPP_

#include <vector>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"  // Should be included before any other Chaste headers.

class ParameterPointData
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
        archive & mQoIs;
        archive & mErrorOccurred;
        archive & mErrorEstimates;
    }

    /** The quantities of interest (QoIs) at this point */
    std::vector<double> mQoIs;

    /** Estimates for errors in the quantities of interest (QoIs) at this point */
    std::vector<double> mErrorEstimates;

    /** Whether an error occurred whilst we were evaluating the QoIs.*/
    bool mErrorOccurred;

    /** Private constructor for archiving only. */
    ParameterPointData(){};

public:
    /**
     * Constructor
     *
     * @param rQoIs  The quantities of interest that were evaluated at this point.
     * @param error  Whether an error occurred during their evaluation.
     */
    ParameterPointData(const std::vector<double>& rQoIs, bool error)
        : mQoIs(rQoIs),
          mErrorOccurred(error)
    {

    }

    /**
     * @return The Quantities of Interest that were evaluated here.
     */
    std::vector<double> GetQoIs() const
    {
        return mQoIs;
    }

    /**
     * A similar method to GetQoIs that returns a const'ed reference
     * to lower overhead of copying data around.
     *
     * @return The Quantities of Interest that were evaluated here.
     */
    const std::vector<double>& rGetQoIs() const
    {
        return mQoIs;
    }

    /**
     * @return The error estimates in the Quantities of Interest that were evaluated here.
     */
    const std::vector<double>& rGetQoIErrorEstimates() const
    {
        if (mErrorEstimates.size()==0u)
        {
            EXCEPTION("Error estimates have not been set on this parameter data point.");
        }
        return mErrorEstimates;
    }

    /**
     * Set the error estimates associated with the QoIs at this point
     * @param errorEstimates  The error estimates.
     */
    void SetErrorEstimates(const std::vector<double>& rErrorEstimates)
    {
        assert(rErrorEstimates.size()==mQoIs.size());
        mErrorEstimates = rErrorEstimates;
    }

    /**
     * @return Whether there are any error estimates for this data point.
     */
    bool HasErrorEstimates()
    {
        return (mErrorEstimates.size() > 0u);
    }

    /**
     * @return  Whether an error occurred during the QoI evaluation.
     */
    bool DidErrorOccur() const
    {
        return mErrorOccurred;
    }

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ParameterPointData)

#endif // PARAMETERPOINTDATA_HPP_
