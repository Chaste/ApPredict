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

#ifndef DOSECALCULATOR_HPP_
#define DOSECALCULATOR_HPP_

#include <vector>

/**
 * Helper class to decide on drug concentrations to test, given command line arguments
 *
 * NB: Now works in 'uM' rather than 'nM'
 */
class DoseCalculator
{
private:
    /** Whether to use pre-specified concentrations or not */
    bool mUseSpecifiedConcs;

    /** Whether to space out test concentrations on a log scale */
    bool mLogScale;

    /**
	 * Either the number of sub-divisions between low and high concentrations
	 * OR
	 * the number of sub-divisions between specified concentrations
	 */
    unsigned mNumSubConcentrations;

    /** Top of a specified concentration range */
    double mTopDose;

    /** Bottom of a specified concentration range */
    double mBottomDose;

    /**
     * The specified concentrations to return
     */
    std::vector<double> mConcentrations;

    /**
	 * Generate a set of concentrations between two values (lowest is never included).
	 *
	 * @param low lowest concentration
	 * @param high highest concentration
	 * @param includeTopDose  whether to add the top concentration to the vector which is returned
	 *
	 * @return a vector of mNumSubConcentrations equally/log spaced concentrations between low and high tests
	 */
    std::vector<double> GetEquallySpacedBetween(double low, double high, bool includeTopDose);

public:
    /**
	 * Constructor for precise test concentrations, to match e.g. an experiment
	 *
	 * @param rPreciseDoses  The exact test concentrations to cover.
	 */
    DoseCalculator(const std::vector<double>& rPreciseDoses);

    /**
	 * Constructor for equally spaced concentrations between a max and min
	 *
	 * @param highDose  maximum test concentration
	 * @param lowDose  minimum test concentration (a control may be added)
	 */
    DoseCalculator(double highDose, double lowDose = 0);

    /**
	 * Constructor to take advantage of command line arguments.
	 */
    DoseCalculator();

    /**
	 * Instruct GetConcentrations() to return entries on a log scale or not
	 * @param logScale whether to spread test concentrations on a log scale.
	 */
    void SetLogScale(bool logScale)
    {
        mLogScale = logScale;
    }

    /**
	 * Specify how many sub-intervals you would like to be returned
	 *
	 * @note GetConcentrations() will add in control (0nM, and 1nM for log scales).
	 */
    void SetNumSubdivisions(unsigned numDivisions)
    {
        mNumSubConcentrations = numDivisions;
    }

    /**
	 * Calculates the drug concentrations to test using the settings applied to the class.
	 *
	 * @return Concentrations at which to test drug action
	 */
    std::vector<double> GetConcentrations(void);
};

#endif // DOSECALCULATOR_HPP_
