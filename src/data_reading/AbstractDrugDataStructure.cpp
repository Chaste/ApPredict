/*

Copyright (c) 2005-2021, University of Oxford.
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

#include "AbstractDrugDataStructure.hpp"
#include "Exception.hpp"

template<unsigned NUM_CHANNELS>
unsigned AbstractDrugDataStructure<NUM_CHANNELS>::GetNumDrugs(void)
{
    return mDrugNames.size();
}

template<unsigned NUM_CHANNELS>
std::string AbstractDrugDataStructure<NUM_CHANNELS>::GetDrugName(unsigned drugIndex)
{
    assert(drugIndex < GetNumDrugs());
    return mDrugNames[drugIndex];
}

template<unsigned NUM_CHANNELS>
unsigned AbstractDrugDataStructure<NUM_CHANNELS>::GetDrugIndex(const std::string& rName)
{
    unsigned idx = UINT_MAX;
    for (unsigned i=0; i<mDrugNames.size(); ++i)
    {
        if (mDrugNames[i] == rName)
        {
            idx = i;
            break;
        }
    }
    if (idx==UINT_MAX)
    {
        EXCEPTION("Drug " << rName << " not found.");
    }
    return idx;
}

template<unsigned NUM_CHANNELS>
double AbstractDrugDataStructure<NUM_CHANNELS>::GetIC50Value(unsigned drugIndex, unsigned channelIndex)
{
    assert(drugIndex < GetNumDrugs());
    assert(channelIndex < NUM_CHANNELS);
    double ic50 = mIc50values[drugIndex](channelIndex);

    /**
     * If an IC50 value takes the value
     *  -1 in the data file that means "effect unknown" - we return a negative IC50 (-1) to show this
     *  -2 in the data file means "known to have no effect" - we return (-2) so that we know to apply no channel block.
     */

    return ic50;
}

template<unsigned NUM_CHANNELS>
double AbstractDrugDataStructure<NUM_CHANNELS>::GetHillCoefficient(unsigned drugIndex, unsigned channelIndex)
{
    assert(drugIndex < GetNumDrugs());
    assert(channelIndex < NUM_CHANNELS);
    double hill = mHillCoefficients[drugIndex](channelIndex);
    if (hill < 0)
    {   // Should default to 1
        hill = 1.0;
    }
    return hill;
}

template<unsigned NUM_CHANNELS>
double AbstractDrugDataStructure<NUM_CHANNELS>::GetSaturationLevel(unsigned drugIndex, unsigned channelIndex)
{
    assert(drugIndex < GetNumDrugs());
    assert(channelIndex < NUM_CHANNELS);
    double sat = mSaturationLevels[drugIndex](channelIndex);

    if (sat < 0)
    {   // Should default to 0%
        sat = 0.0;
    }
    return sat;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractDrugDataStructure<3u>;
template class AbstractDrugDataStructure<4u>;
template class AbstractDrugDataStructure<5u>;
template class AbstractDrugDataStructure<7u>;


