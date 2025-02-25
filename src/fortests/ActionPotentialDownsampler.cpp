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

#include <cmath>
#include <cassert>

#include "ActionPotentialDownsampler.hpp"
#include "OutputFileHandler.hpp"
#include "CommandLineArguments.hpp"
#include "Exception.hpp"

ActionPotentialDownsampler::ActionPotentialDownsampler(const std::string& rFoldername,
                                                       const std::string& rFilename,
                                                       const std::vector<double>& rTimes,
                                                       const std::vector<double>& rVoltages,
                                                       double window,
                                                       double stimulusStart)
{
    OutputFileHandler handler(rFoldername, false); // don't wipe the folder!
    out_stream output_file = handler.OpenOutputFile(rFilename);

    *output_file << "Time(ms)\tMembrane_Voltage(mV)\n";

    double last_voltage_printed = DOUBLE_UNSET;
    double last_time_printed = DOUBLE_UNSET;
    bool printed_last = true; // This is to make sure large jumps print out the step before.

    assert(rTimes.size()>0);
    double start_time_for_this_pace = rTimes[0] + stimulusStart;

    for (unsigned i=0; i<rVoltages.size(); i++)
    {
        if (rTimes[i] - start_time_for_this_pace > window)
        {
        	// Only output the first action potential.
            break;
        }
        if (!CommandLineArguments::Instance()->OptionExists("--no-downsampling"))
        {
            // A new bit of code to downsample the output so flot doesn't timeout loading large datasets
            // We want to plot the first point, but not necessarily the last if we are repolarized.
            if (i>0 && (i<rVoltages.size() || rVoltages[i] < -50)) // if we aren't at the beginning or the end of the trace
            {
            	if (fabs(rVoltages[i] - last_voltage_printed) > 1.0 /*mV*/ )
                {
            		// Don't skip point if the voltage has changed by more than 1 mV
                }
            	else if (rTimes[i] - last_time_printed > 10.0 /*ms*/ )
            	{
            		// Don't skip point if time since last point is over 10ms
            		printed_last = true; // But if this is the reason don't re-print last point too as we do for voltage jumps.
            	}
            	else
            	{
            		// Skip this point
            		printed_last = false;
            		continue;
            	}
            }
            if(!printed_last) // We want the point before printing too, to avoid large interpolations.
            {
            	*output_file << rTimes[i-1] - start_time_for_this_pace << "\t" << rVoltages[i-1] << "\n";
            }
            printed_last = true;
            last_voltage_printed = rVoltages[i];
            last_time_printed = rTimes[i];
        }
        *output_file << rTimes[i] - start_time_for_this_pace << "\t" << rVoltages[i] << "\n";
    }
    output_file->close();
}
