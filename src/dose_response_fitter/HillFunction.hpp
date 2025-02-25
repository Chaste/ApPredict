/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef HILLFUNCTION_HPP_
#define HILLFUNCTION_HPP_

#include <vector>

/**
 * A class that will evaluate an error metric for the match between a Hill function with the specified parameters,
 * and an experimental dose-response curve.
 */
class HillFunction
{
public:
    /**
     * HillFunction constructor
     *
     * @param low  The lowest possible Hill coefficient to assign in the fitting process.
     * @param high  The highest possible Hill coefficient to assign in the fitting process.
     */
	HillFunction(double low, double high);

	/**
	 * To set penalty by which negative and too high hill coefficients will be penalised
	 * @param penalty_value
	 */
	void SetPenalty(double penalty_value);

	/**
	 * Method taking in vector of current parameter values to evaluate error associated with this parameter combination
	 * passed from Nelder-Mead Simplex algorithm after fitting to appropriate hill function (fitting for 1 or 2 parameters,
	 * just IC50 or IC50 and hill coefficient) depending on number of data points input.
	 *
	 * @param rParameters the parameters with which to perform the error function evaluation.
	 * @return Error metric
	 */
	double Evaluate(const std::vector<double>& rParameters);

	/**
	 * Methods to set concentrations as observed in experimental data.
	 * Entries of the vector should have a 1:1 correspondence with inhibitions.
	 *
	 * @param rConcentrations  the doses or concentrations at which responses/inhibitions were recorded.
	 * @param rInhibitions  the % inhibitions that were recorded in the experiment.
	 */
	void SetConcentrationsAndInhibitions(std::vector<double>& rConcentrations, std::vector<double>& rInhibitions);

private:
	// Member variables used in the class

	/** Amount by which a Hill coefficient parameter going out of bounds is penalised - can be set by user.*/
	double mPenalty;
	/** A vector of concentrations at which experimental data points were taken*/
	std::vector<double> mConcentrations;
	/** A vector of inhibitions (expressed as %ages of control current) recorded experimentally.*/
	std::vector<double> mInhibitions;

	/** The maximum value a Hill coefficient should take in the fitting process (defaults to 5 in constructor). */
	double mMinHill;
	/** The minimum value a Hill coefficient should take in the fitting process (defaults to 0 in constructor). */
	double mMaxHill;

};

#endif /* HILLFUNCTION_HPP_ */
