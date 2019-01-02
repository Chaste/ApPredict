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

#ifndef ACTIONPOTENTIALDOWNSAMPLER_HPP_
#define ACTIONPOTENTIALDOWNSAMPLER_HPP_

#include <vector>
#include <string>

/**
 * Class that takes an action potential and downsamples it for good visualization and smaller file size.
 */
class ActionPotentialDownsampler
{
public:
    /**
     * Constructor - this does all the work and outputs a file.
     * Takes an action potential and downsamples it for good visualization and small file size.
     *
     * @param rFoldername  Foldername to write to
     * @param rFilename  Filename to write to
     * @param rTimes  a vector of time points, generally in ms.
     * @param rVoltages  a vector of voltages in mV corresponding to those time points.
     * @param window  the time window to output for (the first 'window' ms will go out to file).
     * @param stimulusStart  the time that the stimulus is applied relative to start of trace (defaults to zero).
     *                       This just allows the upstroke time to be set to zero in the resulting traces if wanted.
     */
    ActionPotentialDownsampler(const std::string& rFoldername,
                               const std::string& rFilename,
                               const std::vector<double>& rTimes,
                               const std::vector<double>& rVoltages,
                               double window,
                               double stimulusStart=0);
};

#endif // ACTIONPOTENTIALDOWNSAMPLER_HPP_
