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

#ifndef LOGISTICDISTRIBUTION_HPP_
#define LOGISTICDISTRIBUTION_HPP_

#include "AbstractDistribution.hpp"

/**
 * This class describes a logistic distribution,
 * it gives access to the probability density function (PDF) and
 * also provides samples from the distribution.
 */
class LogisticDistribution : public AbstractDistribution
{
public:

    /**
     * Constructor
     */
    LogisticDistribution()
	: AbstractDistribution()
	{
	}

    /**
     * Evaluate the probability density function of the Logistic distribution.
     *
     * @param rMu  first (centering) parameter
     * @param rSigma  second (spread) parameter
     * @param rSample  the possible sample value 'X' at which to evaluate the PDF
     * @return PDF(X)
     */
    double EvaluatePdf(const double& rMu, const double& rSigma, const double& rSample)
    {
         return exp(-(rSample - rMu)/rSigma)/(rSigma*(1.0+exp(-(rSample - rMu)/rSigma))*(1.0+exp(-(rSample - rMu)/rSigma)));
    }

private:

	/**
	 * Do the inverse CDF on a given distribution
	 *
	 * @param rMu  parameter 1 of the Logistic
	 * @param rSigma  parameter 2 of the Logistic
	 * @param rP  probability in [0, 1].
	 *
	 * @return the value of X for which CDF(X)=p.
	 */
	double GetSingleSample(const double& rMu, const double& rSigma, const double& rP)
	{
		return rMu + rSigma*log((rP)/(1-rP));
	}
};

#endif // LOGISTICSAMPLER_HPP_
