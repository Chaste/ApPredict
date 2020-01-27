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

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>

#include <boost/scoped_array.hpp>

#include "RunHillFunctionMinimization.hpp"
#include "NelderMeadMinimizer.hpp"

//Concentrations and Inhibitions set inside class constructor
RunHillFunctionMinimization::RunHillFunctionMinimization(const std::vector<double>& rConcentrations,
                                                         const std::vector<double>& rInhibitions,
                                                         unsigned numParameters,
                                                         bool roundValues)
: mConcentrations(rConcentrations),
  mInhibitions(rInhibitions),
  mNumParameters(numParameters),
  mRoundValues(roundValues),
  mTotalNumEvaluations(0u),
  mMinHill(0),
  mMaxHill(5) // 5 is widely considered a sensible limit by GSK and AZ.
{
}

void RunHillFunctionMinimization::SetHillLimits(double low, double high)
{
    mMinHill = low;
    mMaxHill = high;
}

//Method run to prepare the input data - calculating starting parameters for algorithm and selecting
//appropriate hill function to be fitted depending on number of input data points. Parameters are then
//returned here.
std::vector<double> RunHillFunctionMinimization::Run()
{
    mTotalNumEvaluations = 0;
    std::vector<double> parameters;
    if (mNumParameters==2)
    {
        parameters = RunFitForNParams(1);
        assert(parameters.size()==1);
        // Now add an initial guess for Hill coefficient.
        parameters.push_back(1.0);
    }
    parameters = RunFitForNParams(mNumParameters, parameters);
    std::cout << "Minimization complete: total number of function evaluations = " << mTotalNumEvaluations << std::endl << std::flush;

    // A new section to limit the values we get back to sensible limits
    double capping_value = 1e6; // uM
    if (mRoundValues && parameters[0] >= capping_value)
    {
        std::cout << "IC50 that was fitted = " << parameters[0] << "uM, this is outside measurable range.\n\n"
                "So we are capping the fitted value to " << capping_value << " uM.\n(Even at 100uM this corresponds to only 0.01% block)\n"
                "and setting the corresponding Hill coefficient to 1.\n\n";
        // This 1e6 IC50 value corresponds to a pIC50 of 0.
        // At 100uM (a massive concentration) there is still less than 0.01% block at this IC50 (if Hill = 1).
        parameters[0] = capping_value;
        if (mNumParameters==2)
        {   // We may as well set the Hill coefficient to a sensible value too.
            parameters[1] = 1.0;
        }
    }

    return parameters;
}

std::vector<double> RunHillFunctionMinimization::RunFitForNParams(unsigned numParamsToFit, std::vector<double> initialGuess)
{
    // Number of input data points calculated
    unsigned numberofpoints = mConcentrations.size();

    // Number of parameters to be fitted is changed to 1 if all input concentrations are the same
    unsigned repeatedconcs = 0;
    for (unsigned i = 1; i<numberofpoints;i++)
    {
        if (mConcentrations[i]==mConcentrations[i-1])
        {
            repeatedconcs++;
        }
    }
    if (repeatedconcs == numberofpoints-1)
    {
        numParamsToFit = 1;
    }

    // Make a hill function
    HillFunction* p_hill_function = new HillFunction(mMinHill, mMaxHill);

    // Initial guess at parameters for algorithm calculated -
    // concentration at which inhibition recorded
    // is closest to 50% is chosen and initial hill coefficient (where appropriate) is 1
    double ydiffprev = 10000;
    int itrack = -1;

    boost::scoped_array<double> ydiff(new double[numberofpoints]);
    for (unsigned i = 0; i<numberofpoints;i++)
    {
        ydiff[i] = fabs(50 - mInhibitions[i]);
        if (ydiff[i]<ydiffprev)
        {
            itrack = i;

            ydiffprev = ydiff[i];
        }
    }

    // Making sure that the condition imposed to begin the identification process above actually works
    assert(itrack!=-1);

    // Vector of initial guess parameters initialised
    std::vector<double> parameters;
    if (numParamsToFit > 1 && numberofpoints >1)
    {
        if (initialGuess.size()==0)
        {
            parameters.push_back(mConcentrations[itrack]);
            parameters.push_back(1.0);
        }
        else
        {
            parameters = initialGuess;
        }
    }
    else
    {
        if (initialGuess.size()==0)
        {
            parameters.push_back(mConcentrations[itrack]);
        }
        else
        {
            parameters.push_back(initialGuess[0]);
        }
    }

    // Concentrations and Inhibition data, with initial guess at parameters are passed, with the appropriate
    //hill function used to be begin the minimisation.

    //    for (unsigned i=0; i< parameters.size(); i++)
    //    {
    //        std::cout << "InitialCondition[" << i << "] = " << parameters[i] <<"\n";
    //    }

    p_hill_function->SetConcentrationsAndInhibitions(mConcentrations, mInhibitions);
    NelderMeadMinimizer nelder_mead(parameters, p_hill_function);
    //    nelder_mead.SetDisplayIterations(true);
    nelder_mead.Minimize();
    mTotalNumEvaluations += nelder_mead.GetNumEvaluations();

    // deleting from memory
    delete p_hill_function;

    return parameters;
}
