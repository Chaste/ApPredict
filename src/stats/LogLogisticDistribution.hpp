/*

Copyright (c) 2005-2024, University of Oxford.
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

#ifndef LOGLOGISTICDISTRIBUTION_HPP_
#define LOGLOGISTICDISTRIBUTION_HPP_

#include "AbstractDistribution.hpp"

/**
 * This class gives a sample from a log-logistic distribution.
 *
 * NB MatLab parameterises this distribution with the mu and sigma
 * that it would have if we fitted the logistic, i.e.
 * ln(X)~Logistic(mu, sigma)
 *
 * But here we parameterise it with its own parameters
 * X~LogLogistic(alpha,beta)
 * [as in the wikipedia entry as of 14th Oct 2012!]
 *
 * Their relationship is:
 * alpha = exp(mu)
 * beta = 1/sigma
 */
class LogLogisticDistribution : public AbstractDistribution
{
public:

    /**
     * Constructor
     */
    LogLogisticDistribution()
	: AbstractDistribution()
	{
	}

    /**
     * Evaluate the Probability Density Function (PDF) of the LogLogistic.
     *
     * @param rAlpha  The first (centering) parameter of the LogLogistic
     * @param rBeta  The second (spread) parameter of the LogLogistic
     * @param rSample a possible sample 'X'.
     * @return the value of PDF(X)
     */
    double EvaluatePdf(const double& rAlpha, const double& rBeta, const double& rSample)
    {
        double temp = 1.0 + pow(rSample/rAlpha,rBeta);
        return ((rBeta/rAlpha)*pow(rSample/rAlpha,rBeta-1.0))/(temp*temp);
    }

private:
	/**
	 * Get a single sample from the Log-Logistic distribution, with p in [0, 1]
	 *
	 * @param rAlpha - the first, centering, parameter;
	 * @param rBeta - the second, spread, parameter;
	 * @param rP  Probability
	 * @return the value of X for which CDF(X)=p.
	 */
	double GetSingleSample(const double& rAlpha, const double& rBeta, const double& rP)
	{
		return rAlpha*pow(rP/(1.0-rP),1.0/rBeta);
	}
};

#endif // LOGLOGISTICDISTRIBUTION_HPP_
