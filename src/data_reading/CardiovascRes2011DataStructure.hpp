/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef CARDIOVASCRESDATASTRUCTURE_HPP_
#define CARDIOVASCRESDATASTRUCTURE_HPP_

#include "AbstractDrugDataStructure.hpp"

/**
 * Helper class to read in data
 */
class CardiovascRes2011DataStructure: public AbstractDrugDataStructure<3>
{
  protected:
    virtual void LoadALine(std::stringstream& rLine)
    {
        int redfern_category;
        std::string name;
        std::string in_redfern_figs;
        c_vector<double, 3> ic50s;
        ic50s.clear();
        c_vector<double, 2> doses;
        double grandi_measure;
        doses.clear();

        rLine >> name;
        rLine >> redfern_category;
        // INa
        rLine >> ic50s(0);
        // ICaL
        rLine >> ic50s(1);
        // IKr
        rLine >> ic50s(2);
        rLine >> doses(0); // small EFTPC max
        rLine >> doses(1); // large EFTPC max
        rLine >> grandi_measure;

        this->mDrugNames.push_back(name);
        mRedfernCategory.push_back(redfern_category);
        this->mIc50values.push_back(ic50s);
        mClinicalDoseRange.push_back(doses);
        mGrandiMeasure.push_back(grandi_measure);
    }

    std::vector<int> mRedfernCategory;
    std::vector<c_vector<double,2> > mClinicalDoseRange;
    std::vector<double> mGrandiMeasure;

  public:

    CardiovascRes2011DataStructure(std::string fileName)
      : AbstractDrugDataStructure<3>()
    {
        LoadDataFromFile(fileName);
    };

    CardiovascRes2011DataStructure(FileFinder& rFileFinder)
      : AbstractDrugDataStructure<3>()
    {
        LoadDataFromFile(rFileFinder.GetAbsolutePath());
    };

    bool HasRedfernCategory(unsigned drugIndex)
    {
        bool has_redfern = true;
        assert(drugIndex < GetNumDrugs());
        if (mRedfernCategory[drugIndex] < 0)
        {
            has_redfern = false;
        }
        return has_redfern;
    }

    unsigned GetRedfernCategory(unsigned drugIndex)
    {
        assert(drugIndex < GetNumDrugs());
        if (!HasRedfernCategory(drugIndex))
        {   // Not covered because they've all got one at the moment (but might not in the future).
            EXCEPTION("This drug has no Redfern (2003) TdP Liability Category.");
        }
        if (mRedfernCategory[drugIndex]<=0)
        {
            EXCEPTION("Drug " << this->mDrugNames[drugIndex] << " has Redfern category " << mRedfernCategory[drugIndex]);
        }
        assert(mRedfernCategory[drugIndex]>0);
        assert(mRedfernCategory[drugIndex]<6);
        return (unsigned)(mRedfernCategory[drugIndex]);
    }

    double GetGrandiMeasure(unsigned drugIndex)
    {
        assert(drugIndex < GetNumDrugs());
        if (mGrandiMeasure[drugIndex] < -998)
        {
            EXCEPTION("No data available on Grandi measure for " << this->mDrugNames[drugIndex]);
        }
        return mGrandiMeasure[drugIndex];
    }

    double GetClinicalDoseRange(unsigned drugIndex, unsigned lowOrHigh)
    {
        assert(lowOrHigh==0 || lowOrHigh==1);
        assert(drugIndex < GetNumDrugs());
        if (mClinicalDoseRange[drugIndex](0) < 0)
        {
            EXCEPTION("No data available on clinical dose for " << this->mDrugNames[drugIndex]);
        }
        return mClinicalDoseRange[drugIndex](lowOrHigh);
    }

    bool HasClinicalDoseRange(unsigned drugIndex)
    {
        assert(drugIndex < GetNumDrugs());
        bool result = true;
        if (mClinicalDoseRange[drugIndex](0) < 0)
        {
            result = false;
        }
        return result;
    }

};

#endif // CARDIOVASCRESDATASTRUCTURE_HPP_
