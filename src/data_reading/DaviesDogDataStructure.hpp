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

#ifndef DAVIESDOGDATASTRUCTURE_HPP_
#define DAVIESDOGDATASTRUCTURE_HPP_

#include "AbstractDataStructure.hpp"

class DaviesDogDataStructure : public AbstractDataStructure
{
public:
    DaviesDogDataStructure(const FileFinder& rFileFinder);

    unsigned GetNumDogs();

    double GetIKrParameter(unsigned dog);

    double GetItoParameter(unsigned dog);

    double GetINaParameter(unsigned dog);

    double GetICaLParameter(unsigned dog);

    double GetIK1Parameter(unsigned dog);

    double GetICabParameter(unsigned dog);

    double GetIpCaParameter(unsigned dog);

    double GetINcxParameter(unsigned dog);

    double GetINaKParameter(unsigned dog);

    double GetINaLParameter(unsigned dog);

    double GetICaLTauPowerParameter(unsigned dog);

    double GetItoGateParameter(unsigned dog);

private:

    virtual void LoadALine(std::stringstream& rLine);

    std::vector<unsigned> mDogIndices;
    std::vector<double> mIKrFactors;
    std::vector<double> mItoFactors;
    std::vector<double> mINaFactors;
    std::vector<double> mICaLFactors;
    std::vector<double> mIK1Factors;
    std::vector<double> mICabFactors;
    std::vector<double> mIpCaFactors;
    std::vector<double> mINcxFactors;
    std::vector<double> mINaKFactors;
    std::vector<double> mINaLFactors;
    std::vector<double> mTauPowerFactors;
    std::vector<double> mItoGateFactors;

};

#endif // DAVIESDOGDATASTRUCTURE_HPP_
