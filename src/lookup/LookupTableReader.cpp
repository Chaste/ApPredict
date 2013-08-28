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

#include <algorithm>
#include "Exception.hpp"
#include "LookupTableReader.hpp"

template<unsigned DIM>
void LookupTableReader<DIM>::LoadALine(std::stringstream& rLine)
{
    c_vector<double, DIM> params;
    std::vector<double> quantities;

    for (unsigned i=0; i<mNumParameters; i++)
    {
        double number;
        rLine >> number;
        params[i] = number;
    }
    bool err;
    rLine >> err;
    mDidErrorOccur.push_back(err);
    for (unsigned i=0; i<mNumQoIs; i++)
    {
        double number;
        rLine >> number;
        quantities.push_back(number);
    }

    mParameterPoints.push_back(params);
    mQuantitiesOfInterest.push_back(quantities);
}

template<unsigned DIM>
bool LookupTableReader<DIM>::LoadHeaderLine(std::stringstream& rLine)
{
    rLine >> mNumParameters;
    if (mNumParameters!=DIM)
    {
        EXCEPTION("Dimension mismatch: this LookupTableReader<" << DIM <<
                  "> cannot read files of parameter dimension " << mNumParameters << ".");
    }
    rLine >> mNumQoIs;

    for (unsigned i=0; i<mNumParameters; i++)
    {
        std::string param_name;
        rLine >> param_name;
        mParameterNames.push_back(param_name);
    }
    for (unsigned i=0; i<mNumQoIs; i++)
    {
        int qoi;
        rLine >> qoi;
        mQuantitiesToRecord.push_back((QuantityOfInterest)(qoi));
    }
    unsigned num_error_estimates = 0u;
    rLine >> num_error_estimates;
    std::vector<double> errors_these_params;
    errors_these_params.resize(num_error_estimates);
    for (unsigned i=0; i<num_error_estimates; i++)
    {
        rLine >> errors_these_params[i];
    }
    mErrorEstimates.push_back(errors_these_params);

    return true;
}

template<unsigned DIM>
LookupTableReader<DIM>::LookupTableReader(const std::string& rFileName,const std::string& rOutputFolder)
 : AbstractDataStructure()
{
    FileFinder finder(rOutputFolder + "/" + rFileName + ".dat", RelativeTo::ChasteTestOutput);

    LoadDataFromFile(finder.GetAbsolutePath());

    assert(mQuantitiesOfInterest.size()==mParameterPoints.size());
}

template<unsigned DIM>
std::vector<c_vector<double, DIM> > LookupTableReader<DIM>::GetParameterPoints()
{
    return mParameterPoints;
}

template<unsigned DIM>
std::vector<std::vector<double> > LookupTableReader<DIM>::GetFunctionValues()
{
    return mQuantitiesOfInterest;
}

template<unsigned DIM>
std::vector<std::vector<double> > LookupTableReader<DIM>::GetErrorEstimates()
{
    return mErrorEstimates;
}

template<unsigned DIM>
std::vector<QuantityOfInterest> LookupTableReader<DIM>::GetListOfQuantitiesOfInterest()
{
    return mQuantitiesToRecord;
}

template<unsigned DIM>
std::vector<bool> LookupTableReader<DIM>::GetErrors()
{
    return mDidErrorOccur;
}

template<unsigned DIM>
std::vector<double> LookupTableReader<DIM>::GetQuantity(QuantityOfInterest quantity)
{
    bool found = (std::find(mQuantitiesToRecord.begin(), mQuantitiesToRecord.end(), quantity)!=mQuantitiesToRecord.end());
    if (!found)
    {
        EXCEPTION("This quantity of interest was not recorded.");
    }

    unsigned idx;
    for (idx=0; idx<mQuantitiesToRecord.size(); idx++)
    {
        if (mQuantitiesToRecord[idx]==quantity)
        {
            break;
        }
    }

    std::vector<double> qoi;
    std::vector<std::vector<double> >::iterator iter;
    for (iter = mQuantitiesOfInterest.begin();
         iter != mQuantitiesOfInterest.end();
         ++iter)
    {
        qoi.push_back((*iter)[idx]);
    }
    return qoi;
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////
template class LookupTableReader<1>;
template class LookupTableReader<2>;
template class LookupTableReader<3>;
template class LookupTableReader<4>;
template class LookupTableReader<5>;
// Just up to 5D for now, may need to be bigger eventually.



