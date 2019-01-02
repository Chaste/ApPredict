/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef CIPAQNETCALCULATOR_HPP_
#define CIPAQNETCALCULATOR_HPP_

#include <boost/shared_ptr.hpp>

#include "AbstractCvodeCell.hpp"

/**
 * A class to calculate qNet using a pre-steady state paced ORdCiPAv1 model (only).
 * Note that drug block should already have been applied, and the model paced to steady state.
 *
 * Calculation of 'qNet' as per:
 * "Optimization of an In silico Cardiac Cell Model for Proarrhythmia Risk Assessment"
 * Dutta et al. Frontiers in Physiology. https://doi.org/10.3389/fphys.2017.00616
 *
 * Produces a std::numeric_limits<double>::quiet_NaN() if the AP doesn't repolarise.
 */
class CipaQNetCalculator
{
    boost::shared_ptr<AbstractCvodeCell> mpModel;
    boost::shared_ptr<RegularStimulus> mpStimulus;

public:
    /**
	 * Constructor
	 *
	 * @param mpModel  An ORdCipav1 model (this is checked) in steady state for this drug concentration and parameter set.
	 */
    CipaQNetCalculator(boost::shared_ptr<AbstractCvodeCell> pModel);

    /**
     * Destructor (empty)
     */
    virtual ~CipaQNetCalculator(){};

    /**
     * Run one action potential pace at fine resolution and compute qNet.
     *
     * Retuns NaN if there is no full action potential (repolarisation failure).
     *
     * @return qNet -- the integral of the net outward currents during a complete AP.
     */
    double ComputeQNet();
};

#endif /* CIPAQNETCALCULATOR_HPP_ */
