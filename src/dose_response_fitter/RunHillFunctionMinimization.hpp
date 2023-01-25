/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef RUNHILLFUNCTIONMINIMIZATION_HPP_
#define RUNHILLFUNCTIONMINIMIZATION_HPP_
#include <vector>

/**
 * A helper class to perform dose-response curve fitting.
 */
class RunHillFunctionMinimization
{
public:

	/**
	 * Constructor which accepts input concentration and inhibition data vectors and parameter value vector
	 *
	 * @param rConcentrations  a vector of concentrations (does not have to be ordered).
	 * @param rInhibitions  a vector of percentage inhibitions at this concentration (same indices as the concentrations)
	 * @param numParameters  The number of parameters to fit (1 = IC50 only, 2 = IC50 and Hill)
	 * @param roundValues  Whether to round the maximum IC50s that are returned to sensible values (to avoid negative confusing pIC50s).
	 */
	RunHillFunctionMinimization(const std::vector<double>& rConcentrations,
	                            const std::vector<double>& rInhibitions,
	                            unsigned numParameters=2u,
	                            bool roundValues = true);

	/**
	 * Method which takes in concentration and inhibition data vectors, calculates the starting parameters for the
     * algorithm based on this data, decides which form of the hill function needs to be fitted (fitting for
     * 1 or 2 variables depending on number of input data points and then passes these so that the minimization
     * algorithm can go on to be started
     *
     * @return the fitted parameters.
	 */
	std::vector<double> Run();

	/**
	 * Change the limits on the fitting algorithm, so it doesn't return rubbish values (defaults assigned in constructor).
	 *
	 * @param low  The lowest value that a Hill coefficient should be fitted as.
	 * @param high  The highest value that a Hill coefficient should be fitted as.
	 */
	void SetHillLimits(double low, double high);

private:

	/**
	 * Set up initial guesses, make a Hill function, and attempt a fit.
	 * @param numParamsToFit  The number of parameters to fit on this run.
	 * @param initialGuess  An optional initial guess for the minimisation to use.
	 * @return Minimised parameters
	 */
	std::vector<double> RunFitForNParams(unsigned numParamsToFit, std::vector<double> initialGuess = std::vector<double>());

	/** The concentrations at which experimental data points were recorded */
	std::vector<double> mConcentrations;
	/** The inhibitions recorded in experiment */
	std::vector<double> mInhibitions;

	/** The number of parameters to fit (1 = IC50 only, 2 = IC50 and Hill).*/
	unsigned mNumParameters;

	/** Whether to cap the values that can be returned from the Run() method to an IC50 of 1e6 uM = PIC50 of 0.*/
    bool mRoundValues;
    /** The total number of function evaluations that have been made in the fitting process */
	unsigned mTotalNumEvaluations;

    /** The maximum value a Hill coefficient should take in the fitting process (defaults to 5 in constructor). */
    double mMinHill;
    /** The minimum value a Hill coefficient should take in the fitting process (defaults to 0 in constructor). */
    double mMaxHill;

};


#endif /* RUNHILLFUNCTIONMINIMIZATION_HPP_ */
