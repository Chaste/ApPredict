/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef ABSTRACTDISTRIBUTION_HPP_
#define ABSTRACTDISTRIBUTION_HPP_

#include "RandomNumberGenerator.hpp"

/**
 * This class describes all the features that a 'distribution' class should provide,
 * and contains some shared functionality, for example for taking the mean of a number of samples from
 * a concrete distribution.
 */
class AbstractDistribution
{
public:
    /**
     * Constructor (empty)
     */
    AbstractDistribution(){};

    /**
     * Destructor (empty)
     */
	virtual ~AbstractDistribution(){};

	/**
	 * Take a sample at random from a distribution
	 *
	 * @param rParam1 The first parameter of the distribution
	 * @param rParam2 The second parameter of the distribution
	 * @param num_experiments The number of individual samples to take
	 * @return the sample  The mean result of num_experiments samples.
	 */
	double GetSample(const double& rParam1, const double& rParam2, unsigned num_experiments = 1 )
	{
		double return_value = 0.0;
		for (unsigned i=0; i<num_experiments; i++)
		{
			double p = RandomNumberGenerator::Instance()->ranf();
			return_value += GetSingleSample(rParam1, rParam2, p);
		}
		return return_value / ((double)(num_experiments));
	}

	/**
	 * Evaluate the probability density function (PDF) of the distribution at a given point.
	 *
	 * Each sub-class must provide this method.
	 *
	 * @param rParam1  The first parameter of the distribution
	 * @param rParam2  The second parameter of the distribution
	 * @param rSample  The point at which to evaluate the PDF
	 * @return  the PDF evaluated at this point.
	 */
	virtual double EvaluatePdf(const double& rParam1, const double& rParam2, const double& rSample) = 0;

protected:

private:
	/**
	 * This is the inverse of the cumulative density function (CDF),
	 * given a probability p in the range 0<=p<=1, we return the sample
	 * X for which CDF(X) = p. This method must be provided by all subclasses.
	 *
	 * @param rParam1  The first parameter of the distribution
	 * @param rParam2  The second parameter of the distribution
	 * @param rP  the probability in 0<=p<=1
	 * @return  the inverse CDF (value of X for which the CDF(X) = p).
	 */
	virtual double GetSingleSample(const double& rParam1, const double& rParam2, const double& rP) = 0;
};

#endif // ABSTRACTDISTRIBUTION_HPP_
