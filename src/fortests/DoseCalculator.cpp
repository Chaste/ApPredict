/*

Copyright (c) 2005-2014, University of Oxford.
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
#include <algorithm>
#include <cassert>

#include "Exception.hpp"
#include "DoseCalculator.hpp"
#include "CommandLineArguments.hpp"
#include <iostream>

/*
 * CONSTRUCTORS
 */

DoseCalculator::DoseCalculator(const std::vector<double>& rPreciseDoses)
	: mUseSpecifiedConcs(true),
	  mLogScale(false),
	  mNumSubConcentrations(0),
	  mConcentrations(rPreciseDoses)
{
}

DoseCalculator::DoseCalculator(double highDose, double lowDose)
	: mUseSpecifiedConcs(false),
	  mLogScale(false),
	  mNumSubConcentrations(9)
{
	mConcentrations.clear();
	mBottomDose = lowDose;
	mTopDose = highDose;

    // Sanity check
    if (mBottomDose > mTopDose)
    {
        EXCEPTION("Bottom test concentration cannot be larger than top test concentration.");
    }
    assert(mBottomDose>=0);
    assert(mTopDose>0);
}

DoseCalculator::DoseCalculator()
{
    CommandLineArguments* p_args = CommandLineArguments::Instance();
	// First decide whether this is a specified range or precise values
	mUseSpecifiedConcs = p_args->OptionExists("--plasma-concs");
	if (p_args->OptionExists("--plasma-conc-high"))
	{
		mTopDose = p_args->GetDoubleCorrespondingToOption("--plasma-conc-high");
	}
	else
	{
		if (!mUseSpecifiedConcs)
		{
			EXCEPTION("Argument \"--plasma-conc-high <concentration in uM>\" or \"--plasma-concs <concentrations in uM>\" is required");
		}
		else
		{
			mConcentrations = p_args->GetDoublesCorrespondingToOption("--plasma-concs");
		}
	}

	// Get a low value if it exists, zero if not.
	if (p_args->OptionExists("--plasma-conc-low"))
	{
		mBottomDose = p_args->GetDoubleCorrespondingToOption("--plasma-conc-low");
	}
	else
	{
		mBottomDose = 0;
	}

	// Put on a logscale if asked for
	if (p_args->OptionExists("--plasma-conc-logscale"))
	{
		mLogScale = true;
	}
	else
	{
		mLogScale = false;
	}

	// Number of subdivisions to add to the range / precise values...
    if (p_args->OptionExists("--plasma-conc-count"))
    {
    	mNumSubConcentrations = p_args->GetUnsignedCorrespondingToOption("--plasma-conc-count");
    }
    else
    {
    	if (mUseSpecifiedConcs) // Sensible defaults...
    	{
    		mNumSubConcentrations = 0;
    	}
    	else
    	{
    		mNumSubConcentrations = 9;
    	}
    }
}

//
// GENERAL METHODS
//

std::vector<double> DoseCalculator::GetEquallySpacedBetween(double low, double high, bool includeTopDose)
{
	assert(high > low);
	std::vector<double> concs;

	if (mNumSubConcentrations>0)
	{
		for (unsigned i=1 ; i<mNumSubConcentrations + 1; i++)
		{
			if (mLogScale)
			{
				assert(low>0);
				concs.push_back(pow(10,log10(low)+((double)(i)/(double)(mNumSubConcentrations+1))*(log10(high/low))));
			}
			else
			{
				concs.push_back(low + ((double)(i)/(double)(mNumSubConcentrations+1))*(high - low));
			}
		}
	}

	if (includeTopDose)
	{
		concs.push_back(high);
	}

	return concs;
}

std::vector<double> DoseCalculator::GetConcentrations(void)
{
	std::vector<double> concs;
    const double control_for_logscale = 1e-3;

	if (mUseSpecifiedConcs)
	{   // The exact concentrations to test have been input
		concs = mConcentrations;
		std::sort(concs.begin(), concs.end());

		// Make sure there is a control case, if not add it in.
		if (concs[0] > 1e-12)
		{
			concs.push_back(0);
			std::sort(concs.begin(), concs.end());
		}

		// For all Log scale cases also check we have 1nM (1e-3 uM) included
		if (mLogScale && fabs(concs[1]-control_for_logscale) > 1e-12)
		{
			concs.push_back(control_for_logscale);
			std::sort(concs.begin(), concs.end());
		}

		std::vector<double> concs_to_add;
		for (unsigned i=0; i<concs.size()-1; i++)
		{
			if (!(mLogScale && concs[i]==0))
			{
				std::vector<double> concs_to_add_this_subdivision = GetEquallySpacedBetween(concs[i], concs[i+1], false);
				for (unsigned j=0 ; j<concs_to_add_this_subdivision.size() ;j++)
				{
					concs_to_add.push_back(concs_to_add_this_subdivision[j]);
				}
			}
		}
		for (unsigned i=0 ; i<concs_to_add.size() ;i++)
		{
			concs.push_back(concs_to_add[i]);
		}
	}
	else
	{
		/**
		 * DECIDE WHAT SPREAD OF CONCENTRATIONS TO TEST AT
		 */

		// We will always make sure there is a control trace

		// if bottom dose is not zero we add zero in too
		if (fabs(mBottomDose) > 1e-12)
		{
			concs.push_back(0);
		}

		// Put in the lowest dose (zero if no low dose specified).
		concs.push_back(mBottomDose); // The GetEquallySpacedBetween() call below doesn't add the lowest dose.

		// If we are using a log scale and the bottom dose is greater than 1nM we add 1nM too (1e-3 uM).
		if (mLogScale && fabs(mBottomDose-control_for_logscale) > 1e-12)
		{
			concs.push_back(control_for_logscale);
			if (mBottomDose <= control_for_logscale)
			{
				mBottomDose = control_for_logscale;
			}
//			// This snippet of code puts spaced points from control (or 1nM in log) up to the bottom dose
//			// I've decided this can be done better by setting control as bottom dose (or omitting bottom dose) if anyone wants this function.
//			else
//			{
//				std::vector<double> concs_to_add = GetEquallySpacedBetween(control_for_logscale, mBottomDose, false);
//
//				for (unsigned i=0 ; i<concs_to_add.size() ;i++)
//				{
//					concs.push_back(concs_to_add[i]);
//				}
//			}
		}

		// Add in entries between low dose (not included) and top dose (included)
		std::vector<double> concs_to_add = GetEquallySpacedBetween(mBottomDose, mTopDose, true);
		for (unsigned i=0 ; i<concs_to_add.size() ;i++)
		{
			concs.push_back(concs_to_add[i]);
		}
	}
	std::sort(concs.begin(), concs.end());

//	std::cout << "DOSES CALCULATED:\n";
//	for (unsigned i=0; i<concs.size() ; i++)
//	{
//		std::cout << "concs[" << i << "] = " << concs[i] << "\n" << std::flush;
//	}
	return concs;
}



