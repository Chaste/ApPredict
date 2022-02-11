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

#ifndef SINGLEACTIONPOTENTIALPREDICTION_HPP_
#define SINGLEACTIONPOTENTIALPREDICTION_HPP_

#include "AbstractActionPotentialMethod.hpp"
#include "CipaQNetCalculator.hpp"

/**
 * A class to be used to construct lookup tables for drug action,
 * it applies to a model that has been given a scaling factor and
 * initial conditions, runs to steady state, evaluates action
 * potential properties, and returns these.
 */
class SingleActionPotentialPrediction : public AbstractActionPotentialMethod
{
public:
    /**
     * Constructor
     */
    SingleActionPotentialPrediction(boost::shared_ptr<AbstractCvodeCell> pModel)
            : AbstractActionPotentialMethod(),
              mApd90(DOUBLE_UNSET),
              mApd50(DOUBLE_UNSET),
              mUpstroke(DOUBLE_UNSET),
              mPeak(DOUBLE_UNSET),
              mPeakTime(DOUBLE_UNSET),
              mCaMin(DOUBLE_UNSET),
              mCaMax(DOUBLE_UNSET),
              mpModel(pModel){};

    /**
     * Run the steady state pacing and evaluate action potential markers (APD90, APD50, Peak and upstroke velocity).
     * @return  an action potential trace
     */
    OdeSolution RunSteadyPacingExperiment()
    {
        double printing_timestep = 0.1;
        return AbstractActionPotentialMethod::SteadyStatePacingExperiment(mpModel,
                                                                          mApd90,
                                                                          mApd50,
                                                                          mUpstroke,
                                                                          mPeak,
                                                                          mPeakTime,
                                                                          mCaMax,
                                                                          mCaMin,
                                                                          printing_timestep);
    }

    /**
     * Run the steady state pacing and evaluate action potential markers (APD90, APD50, Peak and upstroke velocity).
     *
     * @param conc  The concentration of drug (for more intelligible error messages).
     * @return  an action potential trace
     */
    OdeSolution RunSteadyPacingExperiment(double conc)
    {
        double printing_timestep = 0.1;
        return AbstractActionPotentialMethod::SteadyStatePacingExperiment(mpModel,
                                                                          mApd90,
                                                                          mApd50,
                                                                          mUpstroke,
                                                                          mPeak,
                                                                          mPeakTime,
                                                                          mCaMax,
                                                                          mCaMin,
                                                                          printing_timestep,
                                                                          conc);
    }

    /**
     * @return The APD90 of the trace (ms).
     */
    double GetApd90()
    {
        CheckItRan();
        return mApd90;
    }

    /**
     * @return The APD50 of the trace (ms).
     */
    double GetApd50()
    {
        CheckItRan();
        return mApd50;
    }

    /**
     * @return The qNet for the trace (ms).
     */
    double CalculateQNet()
    {
        CheckItRan();

        CipaQNetCalculator calculator(mpModel);
        return calculator.ComputeQNet();
    }

    /**
     * @return the velocity of the upstroke (mV/ms).
     */
    double GetUpstrokeVelocity()
    {
        CheckItRan();
        return mUpstroke;
    }

    /**
     * @return The peak voltage of the action potential (mV).
     */
    double GetPeakVoltage()
    {
        CheckItRan();
        return mPeak;
    }

    /**
     * @return The time at which peak voltage occurred during the action potential (mV).
     */
    double GetTimeOfPeakVoltage()
    {
        CheckItRan();
        return mPeakTime;
    }

    /**
     * @return The diastolic (minimum) of Calcium transient (mM)
     */
    double GetCaMin()
    {
        CheckItRan();
        return mCaMin;
    }

    /**
     * @return The sysstolic (maximum) of Calcium transient (mM)
     */
    double GetCaMax()
    {
        CheckItRan();
        return mCaMax;
    }

    /**
     * @param pModel a cell model
     *
     * @return the threshold at which we think a voltage signal is a real
     * excited action potential, rather than simply a stimulus current and decay
     * so we give the error code NoActionPotential_1 appropriately (rather than
     * really small APDs).     *
     */
    double DetectVoltageThresholdForActionPotential()
    {
        OdeSolution baseline_solution = RunSteadyPacingExperiment();
        std::vector<double> baseline_voltages = baseline_solution.GetAnyVariable("membrane_voltage");
        double max_baseline_voltage = *(std::max_element(baseline_voltages.begin(), baseline_voltages.end()));
        double min_baseline_voltage = *(std::min_element(baseline_voltages.begin(), baseline_voltages.end()));

        std::string fast_sodium_name = "membrane_fast_sodium_current_conductance";
        std::string l_type_cal_name = "membrane_L_type_calcium_current_conductance";
        if (mpModel->HasParameter(fast_sodium_name + "_scaling_factor"))
        {
            fast_sodium_name += "_scaling_factor";
        }
        if (mpModel->HasParameter(l_type_cal_name + "_scaling_factor"))
        {
            l_type_cal_name += "_scaling_factor";
        }

        // We switch off the sodium and calcium currents and see how high the stimulus makes the voltage go.
        if (mpModel->HasParameter(fast_sodium_name))
        {
            // Remember original conductances and set them to zero.
            const double original_na_conductance = mpModel->GetParameter(fast_sodium_name);
            double original_ca_conductance = 0;
            mpModel->SetParameter(fast_sodium_name, 0u);
            if (mpModel->HasParameter(l_type_cal_name))
            {
                original_ca_conductance = mpModel->GetParameter(l_type_cal_name);
                mpModel->SetParameter(l_type_cal_name, 0u);
            }

            // Remember state variables
            N_Vector steady_full_conductance_state_vars = mpModel->GetStateVariables();

            OdeSolution solution = RunSteadyPacingExperiment();

            // Put conductances back where they were!
            mpModel->SetParameter(fast_sodium_name, original_na_conductance);
            if (mpModel->HasParameter(l_type_cal_name))
            {
                mpModel->SetParameter(l_type_cal_name, original_ca_conductance);
            }

            mpModel->SetStateVariables(steady_full_conductance_state_vars);
            DeleteVector(steady_full_conductance_state_vars);

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

private:
    /**
     * Check whether the run has completed and was unsuccessful,
     * or has not been run yet,
     * throw an exception if either of these is the case.
     */
    void CheckItRan()
    {
        if (!(this->mSuccessful))
        {
            EXCEPTION("We have not run to steady state yet, or the marker evaluation failed.");
        }
    }

    double mApd90;
    double mApd50;
    double mUpstroke;
    double mPeak;
    double mPeakTime;
    double mCaMin;
    double mCaMax;
    boost::shared_ptr<AbstractCvodeCell> mpModel;
};

#endif // SINGLEACTIONPOTENTIALPREDICTION_HPP_
