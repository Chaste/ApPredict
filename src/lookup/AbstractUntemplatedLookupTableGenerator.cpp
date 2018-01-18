/*

Copyright (c) 2005-2018, University of Oxford.
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

#include "AbstractUntemplatedLookupTableGenerator.hpp"

double AbstractUntemplatedLookupTableGenerator::DetectVoltageThresholdForActionPotential(
    boost::shared_ptr<AbstractCvodeCell> pModel)
{
    SingleActionPotentialPrediction ap_runner(pModel);
    ap_runner.SuppressOutput();
    ap_runner.SetMaxNumPaces(100u);

    OdeSolution baseline_solution = ap_runner.RunSteadyPacingExperiment();
    std::vector<double> baseline_voltages = baseline_solution.GetAnyVariable("membrane_voltage");
    double max_baseline_voltage = *(std::max_element(baseline_voltages.begin(), baseline_voltages.end()));
    double min_baseline_voltage = *(std::min_element(baseline_voltages.begin(), baseline_voltages.end()));

    // We switch off the sodium current and see how high the stimulus makes the
    // voltage go.
    if (pModel->HasParameter("membrane_fast_sodium_current_conductance"))
    {
        const double original_na_conductance = pModel->GetParameter("membrane_fast_sodium_current_conductance");
        pModel->SetParameter("membrane_fast_sodium_current_conductance", 0u);

        OdeSolution solution = ap_runner.RunSteadyPacingExperiment();

        // Put it back where it was! The calling method will reset state variables.
        pModel->SetParameter("membrane_fast_sodium_current_conductance",
                             original_na_conductance);

        std::vector<double> voltages = solution.GetAnyVariable("membrane_voltage");
        double max_voltage = *(std::max_element(voltages.begin(), voltages.end()));
        double min_voltage = *(std::min_element(voltages.begin(), voltages.end()));

        // Go 25% over the depolarization jump at gNa=0 as a threshold for 'this
        // really is an AP'. This should be sensible for all models that fail to depolarise with gNa=0.
        const double proposed_threshold = min_voltage + 1.25 * (max_voltage - min_voltage);

        // BUT some models with a big stimulus fire off almost fully with gNa=0 anyway (e.g. ten Tusscher 2006)
        // and so we need to prevent this proposed threshold being above (or near) the usual peak voltage
        // (if we set threshold above we'd always get NoAP1 error codes, even when there are APs!)
        const double two_thirds_of_full_AP = min_baseline_voltage + 0.666 * (max_baseline_voltage - min_baseline_voltage);
        if (proposed_threshold <= two_thirds_of_full_AP)
        {
            return proposed_threshold;
        }
    }

    // Otherwise we give a sensible default of 1/3 of the way up an AP.
    return min_baseline_voltage + 0.333 * (max_baseline_voltage - min_baseline_voltage); // mV
}
