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

// Chaste includes
#include "CellProperties.hpp"
#include "Warnings.hpp"

// ApPredict includes
#include "CipaQNetCalculator.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"

CipaQNetCalculator::CipaQNetCalculator(boost::shared_ptr<AbstractCvodeCell> pModel)
        : mpModel(pModel)
{
    if (!boost::dynamic_pointer_cast<Cellohara_rudy_cipa_v1_2017FromCellMLCvode>(mpModel))
    {
        EXCEPTION("Model used in constructor of CipaQNetCalculator needs to be an Ohara-Rudy-CiPA-v1 model.");
    }

    mpStimulus = boost::static_pointer_cast<RegularStimulus>(mpModel->GetStimulusFunction());
    if (mpStimulus->GetPeriod() != 2000)
    {
        WARNING("qNet should be calculated at 0.5Hz (pacing cycle length of 2000ms), your stimulus is set to " << mpStimulus->GetPeriod() << "ms.");
    }
}

double CipaQNetCalculator::ComputeQNet()
{
    const double maximum_time_step = mpStimulus->GetDuration();
    const double sampling_interval = 0.01; //ms
    OdeSolution detailed_solution = mpModel->Solve(0, mpStimulus->GetPeriod(), maximum_time_step, sampling_interval);

    std::vector<double> times = detailed_solution.rGetTimes();
    std::vector<double> voltages = detailed_solution.GetAnyVariable("membrane_voltage");

    // Analyse trace to check there is actually an AP.
    CellProperties cell_properties(voltages, times, 0 /* mV threshold */);
    try
    {
        cell_properties.GetTimesAtMaxUpstrokeVelocity();
        cell_properties.GetAllActionPotentialDurations(30);
        cell_properties.GetAllActionPotentialDurations(90);
    }
    catch (Exception& e)
    {
        if (e.GetShortMessage() == "AP did not occur, never descended past threshold voltage.")
        {
            std::cout << "Repolarisation Failure, logging NaN for the qNet calculation" << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            throw(e);
        }
    }

    detailed_solution.CalculateDerivedQuantitiesAndParameters(mpModel.get());

    // Inet = ICaL + INaL + IKr + IKs + IK1 + Ito
    std::vector<double> i_cal = detailed_solution.GetAnyVariable("membrane_L_type_calcium_current");
    std::vector<double> i_nal = detailed_solution.GetAnyVariable("membrane_persistent_sodium_current");
    std::vector<double> i_kr = detailed_solution.GetAnyVariable("membrane_rapid_delayed_rectifier_potassium_current");
    std::vector<double> i_ks = detailed_solution.GetAnyVariable("membrane_slow_delayed_rectifier_potassium_current");
    std::vector<double> i_k1 = detailed_solution.GetAnyVariable("membrane_inward_rectifier_potassium_current");
    std::vector<double> i_to = detailed_solution.GetAnyVariable("membrane_transient_outward_current");
    std::vector<double> i_net(i_cal.size());

    // assume equal spacing and convert from milliseconds to seconds.
    const double timestep_in_seconds = (detailed_solution.rGetTimes()[1] - detailed_solution.rGetTimes()[0]) / 1000.0;

    // Integrate over the whole AP.
    double q_net = 0.0; // Charge
    i_net[0] = i_cal[0] + i_nal[0] + i_kr[0] + i_ks[0] + i_k1[0] + i_to[0]; // Current

    for (unsigned i = 1u; i < i_cal.size(); i++)
    {
        i_net[i] = i_cal[i] + i_nal[i] + i_kr[i] + i_ks[i] + i_k1[i] + i_to[i]; // uA/uF

        // Hack in a trapezium rule for the integral, should be enough accuracy with this time resolution.
        q_net += timestep_in_seconds * 0.5 * (i_net[i] + i_net[i - 1]); // Coulomb is Amps * seconds (rather than milliseconds), so integral in uC/uF.
    }

    // Integrate to get qNet
    //std::cout << "Total qNet = " << q_net << " C/F" << std::endl;

    return q_net;
}
