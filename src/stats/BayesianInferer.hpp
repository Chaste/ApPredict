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

#ifndef BAYESIANINFERER_HPP_
#define BAYESIANINFERER_HPP_

#include "AbstractDistribution.hpp"
#include "DoseResponseParameterTypes.hpp"

/**
 * This class works with logistic and log-logistic distributions.
 * It assumes that we have an underlying distribution with known spread (beta or sigma)
 * and unknown median (mu or alpha).
 *
 * By comparing PDF for lots of possible values for mu with the observed data it
 * builds a posterior distribution for the underlying median of the real distribution.
 *
 * We then sample from this to use as possibilities for 'true' dose-response parameter values.
 */
class BayesianInferer
{
public:
  /**
     * Constructor
     *
     * This sets up the possible range of mu (logistic) or alpha (log-logistic) values
     * that the inferer is going to consider as its uniform prior.
     *
     * @param parameter  Which dose-response curve parameter (dictates
     *                   what underlying distribution) this inferer is going to work with.
     */
  BayesianInferer(DoseResponseParameter parameter);

  /**
     * Destructor - just cleans up memory
     */
  ~BayesianInferer();

  /**
     * Set the observed data that we are going to use for inference.
     *
     * These data should stay alive in memory externally,
     * as this class just takes a pointer to them.
     *
     * @param rData the data points
     */
  void SetObservedData(const std::vector<double> &rData);

  /**
     * Set the spread parameter ('sigma' for logistic or 'beta' for log-logistic)
     * @param spread  the spread parameter to assume the underlying distributions have.
     */
  void SetSpreadOfUnderlyingDistribution(double spread);

  /**
     * Get the spread parameter ('sigma' for logistic or 'beta' for log-logistic)
     * @return  the spread parameter to assume the underlying distributions have.
     */
  double GetSpreadOfUnderlyingDistribution();

  /**
     * Perform the inference calculations
     */
  void PerformInference();

  /**
     * Take a sample from the inferred probability distribution.
     *
     * @return a sample of a possible median value from the inferred underlying distribution (of distributions!).
     */
  double GetSampleMedianValue();

  /**
     * This method simply calls GetSampleMedianValue() a number of times and returns you all the values in a vector
     * (to minimise re-calculation), it has nothing to do with repeated data points!
     *
     * @param numValues  the number of samples to take.
     * @return a vector of a possible median values from the inferred underlying distribution (of distributions!).
     */
  std::vector<double> GetSampleMedianValue(const unsigned numValues);

  /**
     * @return The possible Median values (useful for plotting)
     */
  std::vector<double> GetPossibleMedianValues();

  /**
     * @return the posterior CDF (corresponding to the possible median values, useful for plotting).
     */
  std::vector<double> GetPosteriorCdf();

  /**
     * @return the posterior PDF (corresponding to the possible median values, useful for plotting).
     */
  std::vector<double> GetPosteriorPdf();

private:
  /** The parameter we are doing inference on, just an enumeration defined in DoseResponseParameterTypes.hpp*/
  DoseResponseParameter mParameter;

  /** The spread of the underlying distributions */
  double mSigma;

  /** The distribution we are using */
  AbstractDistribution *mpDistribution;

  /** Whether we are ready to perform inference */
  bool mInferenceReady;

  /** The observed data we're working with */
  const std::vector<double> *mpData;

  /** The range of possible median values that we are considering (set in constructor) */
  std::vector<double> mPossibleMuValues;

  /** The posterior PDF */
  std::vector<double> mPosteriorPdf;

  /** The posterior CDF (this is used to take samples)*/
  std::vector<double> mPosteriorCdf;
};

#endif // BAYESIANINFERER_HPP_
