/*

Copyright (c) 2005-2018, University of Oxford.
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

#ifndef ABSTRACTDRUGDATASTRUCTURE_HPP_
#define ABSTRACTDRUGDATASTRUCTURE_HPP_

#include "AbstractDataStructure.hpp"

/**
 * A class which is designed to read in data on drugs from a text file,
 * subclassed for particular formats.
 *
 * This is for any text file with a drug on each line.
 */
template<unsigned NUM_CHANNELS>
class AbstractDrugDataStructure : public AbstractDataStructure
{
protected:

    std::vector<std::string> mDrugNames;
    std::vector<c_vector<double,NUM_CHANNELS> > mIc50values;
    std::vector<c_vector<double,NUM_CHANNELS> > mHillCoefficients;
    std::vector<c_vector<double,NUM_CHANNELS> > mSaturationLevels;

public:

    /**
     * Default Constructor (empty)
     */
    AbstractDrugDataStructure()
     : AbstractDataStructure()
    {};

    /**
     * Destructor (empty)
     */
    virtual ~AbstractDrugDataStructure(){};

    /**
     * @return The number of drugs in the data file.
     */
    unsigned GetNumDrugs(void);

    /**
     * @param drugIndex  The index of the drug (row of the data file on which it appears)
     * @return  The name of the drug
     */
    std::string GetDrugName(unsigned drugIndex);

    /**
     * @param rName  The name of the drug
     * @return  The index in the current drug list.
     */
    unsigned GetDrugIndex(const std::string& rName);

    /**
     *
     * Return the IC50 value associated with a particular channel.
     *
     * If the IC50 is not known then -1 is returned.
     *
     * If the IC50 is known and there is no effect -2 is returned.
     *
     * @param drugIndex  index of the drug as listed in rows of the data file.
     * @param channelIndex  index of the channel we are interested in.
     * @return the IC50 value for a certain drug on a certain channel.
     */
    double GetIC50Value(unsigned drugIndex, unsigned channelIndex);

    /**
     * Return the hill coefficient associated with this drug and this channel's dose-reponse curve
     *
     * @param drugIndex  The index of the drug (in drug_data.dat file)
     * @param channelIndex  The index of the channel we are interested in.
     * @return the hill coefficient.
     */
    double GetHillCoefficient(unsigned drugIndex, unsigned channelIndex);

    /**
     * Return the saturation levels associates with this channel's dose-response curve.
     * Defined relative to 100% - no drug effect. So a value of zero is a drug which can
     * fully block a channel (usually the default),
     * and a value of 150% could activate the channel by an extra 50%.
     *
     * @param drugIndex  The index of the drug (in drug_data.dat file)
     * @param channelIndex  The index of the channel we are interested in.
     * @return the saturation level.
     */
    double GetSaturationLevel(unsigned drugIndex, unsigned channelIndex);

};

#endif // ABSTRACTDRUGDATASTRUCTURE_HPP_
