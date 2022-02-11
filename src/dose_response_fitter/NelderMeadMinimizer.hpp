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

#ifndef NELDERMEADMINIMIZER_HPP_
#define NELDERMEADMINIMIZER_HPP_
#include "HillFunction.hpp"
#include <vector>

class NelderMeadMinimizer
{
public:

	/**
	 * Constructor which accepts input parameter value vector
	 * passed by reference and pointer to appropriate hill function to be used
	 */
	NelderMeadMinimizer(std::vector<double>& rParameters, HillFunction* p_hill_function);

	/**
	 * To set maximum number of iterations Nelder-Mead algorithm can take
	 */
	void SetMaxNumIterations(unsigned numIterations)
	{
		mMaxNumIterations = numIterations;
	}

	/**
	 * @return the number of function evaluations that were performed.
	 */
	unsigned GetNumEvaluations(void)
	{
	    return mNumFunctionEvaluations;
	}

	/**
	 * To set tolerance specified for acceptance of parameters yielding minimum error
	 *
	 * @param tolerance  the tolerance of ?? at which to stop iterating.
	 */
	void SetTolerance(double tolerance)
	{
		mTolerance = tolerance;
	}

	/**
	 * To set reflection coefficient used in Nelder-Mead algorithm
	 *
	 * @param reflectionCoefficient  The reflection coefficient to use.
	 */
	void SetReflectionCoeff(double reflectionCoefficient)
	{
		mReflectionCoefficient = reflectionCoefficient;
	}

	/**
	 * To set expansion coefficient used in Nelder-Mead algorithm
	 *
	 * @param expansionCoefficient  The expansion coefficient to use.
	 */
	void SetExpansionCoeff(double expansionCoefficient)
	{
		mExpansionCoefficient = expansionCoefficient;
	}

	/**
	 * To set contraction coefficient used in Nelder-Mead algorithm
	 *
	 * @param contrationCoefficient  The contraction coefficient to use.
	 */
	void SetContractionCoeff(double contrationCoefficient)
	{
		mContractionCoefficient = contrationCoefficient;
	}

	/**
	 * To set shrink coefficient used in Nelder-Mead algorithm
	 *
	 * @param shrinkCoefficient  The shrink coefficient to use.
	 */
	void SetShrinkCoeff(double shrinkCoefficient)
	{
		mShrinkCoefficient = shrinkCoefficient;
	}

	/**
	 * Allows information on iterations and convergence to be output to screen
	 *
	 * @param display Whether to print information on the iterations to std::cout
	 */
	void SetDisplayIterations(bool display=true)
	{
	    mDisplayIterations = display;
	}

	/**
	 * Method in which minimization algorithm using the Nelder-Mead simplex method is implemented
	 * takes in parameters, calculates new parameters based on rules described in algorithm and
	 * and then passes this to hill function to calculate error associated with fitting with this
	 * parameters. This is repeated until the returned error for parameter combination converges.
	 */
	void Minimize();

private:
	//Reference to vector of parameter values to be passed around
	std::vector<double>& mrParameters;

	//Pointer to appropriate hill function (fitting for one or 2 parameters depending on number of input data points)
	HillFunction* mpFunctionToMinimise;

	//member variables used in the methods in class
	unsigned mNumParameters;

	//member variables which can be set by the user. If not set values set as default in constructor are used
	unsigned mMaxNumIterations;
	double mReflectionCoefficient;
	double mTolerance;
	double mContractionCoefficient;
	double mExpansionCoefficient;
	double mShrinkCoefficient;

	bool mDisplayIterations;

	/** Keep track of the number of function evaluations we perform */
	unsigned mNumFunctionEvaluations;

};

#endif /* NELDERMEADMINIMIZER_HPP_ */
