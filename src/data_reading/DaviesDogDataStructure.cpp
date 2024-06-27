/*

Copyright (c) 2005-2024, University of Oxford.
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

#include "DaviesDogDataStructure.hpp"

DaviesDogDataStructure::DaviesDogDataStructure(const FileFinder& rFileFinder)
  : AbstractDataStructure()
{
    LoadDataFromFile(rFileFinder.GetAbsolutePath());

    assert(mDogIndices.size()==20u);
    for (unsigned i=0; i<mDogIndices.size(); i++)
    {
        assert(i==mDogIndices[i]);
    }
};

unsigned DaviesDogDataStructure::GetNumDogs()
{
    return mDogIndices.size();
}

double DaviesDogDataStructure::GetIKrParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mIKrFactors[dog];
}

double DaviesDogDataStructure::GetItoParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mItoFactors[dog];
}

double DaviesDogDataStructure::GetINaParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mINaFactors[dog];
}

double DaviesDogDataStructure::GetICaLParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mICaLFactors[dog];
}


double DaviesDogDataStructure::GetIK1Parameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mIK1Factors[dog];
}

double DaviesDogDataStructure::GetICabParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mICabFactors[dog];
}

double DaviesDogDataStructure::GetIpCaParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mIpCaFactors[dog];
}

double DaviesDogDataStructure::GetINcxParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mINcxFactors[dog];
}

double DaviesDogDataStructure::GetINaKParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mINaKFactors[dog];
}

double DaviesDogDataStructure::GetINaLParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mINaLFactors[dog];
}

double DaviesDogDataStructure::GetICaLTauPowerParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mTauPowerFactors[dog];
}

double DaviesDogDataStructure::GetItoGateParameter(unsigned dog)
{
    assert(dog < mDogIndices.size());
    return mItoGateFactors[dog];
}

void DaviesDogDataStructure::LoadALine(std::stringstream& rLine)
{
    unsigned dog_index;
    rLine >> dog_index;
    mDogIndices.push_back(dog_index);
    double temp;
    rLine >> temp;
    mIKrFactors.push_back(temp);
    rLine >> temp;
    mItoFactors.push_back(temp);
    rLine >> temp;
    mINaFactors.push_back(temp);
    rLine >> temp;
    mICaLFactors.push_back(temp);
    rLine >> temp;
    mIK1Factors.push_back(temp);
    rLine >> temp;
    mICabFactors.push_back(temp);
    rLine >> temp;
    mIpCaFactors.push_back(temp);
    rLine >> temp;
    mINcxFactors.push_back(temp);
    rLine >> temp;
    mINaKFactors.push_back(temp);
    rLine >> temp;
    mINaLFactors.push_back(temp);
    rLine >> temp;
    mTauPowerFactors.push_back(temp);
    rLine >> temp;
    mItoGateFactors.push_back(temp);
};

