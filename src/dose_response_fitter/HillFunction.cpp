/*

Copyright (c) 2005-2022, University of Oxford.
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

#include <cmath>
#include <cassert>
#include <vector>
#include "HillFunction.hpp"

//HillFunction constructor with default variables for initialisation
HillFunction::HillFunction(double low, double high)
    : mPenalty(1e10),
      mMinHill(low),
      mMaxHill(high)
{
}

//To set penalty by which negative and too high hill coefficients will be penalised
void HillFunction::SetPenalty(double penalty_value)
{
    mPenalty = penalty_value;
}

//Methods to set concentrations and inhibitions from initial input data
void HillFunction::SetConcentrationsAndInhibitions(std::vector<double>& rConcentrations, std::vector<double>& rInhibitions)
{
    assert(rConcentrations.size()==rInhibitions.size());
    mConcentrations = rConcentrations;
    mInhibitions = rInhibitions;
}

// Evaluates error after fitting to appropriate hill function with parameters passed from
// Nelder-Mead algorithm used

double HillFunction::Evaluate(const std::vector<double>& rParameters)
{
	//Assigning hill coefficient appropriately depending on whether we are fitting for 1 or 2 parameters

	double hill_coefficient = 1.0;
	if (rParameters.size() > 1)
	{
		hill_coefficient = rParameters[1];
	}

	double cumerror = 0; // Cumulative error

	for (unsigned i=0; i<mConcentrations.size(); i++)
	{
		double expected_inhibition = (100.0/(1.0+(pow((rParameters[0]/mConcentrations[i]),hill_coefficient))));

		double error = expected_inhibition - mInhibitions[i];

		cumerror += error*error;
	}

    // Penalising error for negative ic50 according to specified penalty
    if (rParameters[0] < 0)
    {
        cumerror -= mPenalty*rParameters[0];
    }

    if (rParameters.size() > 1) // If we are fitting Hill coefficient
    {
        // Penalising error for too small a hill coefficient according to specified penalty
        if (hill_coefficient < mMinHill)
        {
            cumerror += mPenalty*(mMinHill - hill_coefficient);
        }

        //Penalising for too high a hill coefficient according to specified penalty
        if(hill_coefficient > mMaxHill)
        {
            cumerror += mPenalty*(hill_coefficient - mMaxHill);
        }
    }

	return (cumerror);
}
