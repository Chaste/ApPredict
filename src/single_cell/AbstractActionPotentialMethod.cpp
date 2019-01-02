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

// Standard headers
#include <cmath>

// Chaste includes
#include "SteadyStateRunner.hpp"

// ApPredict includes
#include "AbstractActionPotentialMethod.hpp"

AbstractActionPotentialMethod::AbstractActionPotentialMethod()
        : mRunYet(false),
          mMaxNumPaces(UNSIGNED_UNSET),
          mErrorMessage(""),
          mErrorCode(0u),
          mNoOneToOneCorrespondenceIsError(false),
          mAlternansIsError(false),
          mActionPotentialThreshold(-50),
          mActionPotentialThresholdSetManually(false),
          mDefaultParametersApd90(DOUBLE_UNSET),
          mDefaultParametersTimeOfVMax(DOUBLE_UNSET),
          mRepeat(false),
          mRepeatNumber(0u),
          mSuppressOutput(false),
          mHertz(1.0), // default to 1 Hz, replaced by suitable command line
          // argument if present.
          mSuccessful(false),
          mPeriodTwoBehaviour(false)
{
    CommandLineArguments* p_args = CommandLineArguments::Instance();

    if (p_args->OptionExists("--pacing-freq"))
    {
        mHertz = p_args->GetDoubleCorrespondingToOption("--pacing-freq");
        if (mHertz <= 0 + DBL_MIN)
        {
            EXCEPTION(
                "The pacing frequency ("
                << mHertz
                << ") set by '--pacing-freq' option must be a positive number.");
        }
    }
    // mSuppressOutput is always set to false in this constructor, so suppress
    // anyway!
    // if (!mSuppressOutput) std::cout << "* Pacing Frequency = " << mHertz << "
    // Hz\n";

    if (p_args->OptionExists("--pacing-max-time"))
    {
        // Max number of paces is equal to maximum pacing time (converted from
        // minutes to seconds) x paces per second.
        unsigned max_num_paces = std::ceil(60.0 * mHertz * p_args->GetDoubleCorrespondingToOption("--pacing-max-time"));
        SetMaxNumPaces(max_num_paces);
        // if (!mSuppressOutput) std::cout << "* Max Number of paces = " <<
        // max_num_paces << " Hz\n";
    }
}

void AbstractActionPotentialMethod::SetMaxNumPaces(unsigned numPaces)
{
    if (numPaces != UNSIGNED_UNSET)
    {
        mMaxNumPaces = numPaces;
    }
}

unsigned AbstractActionPotentialMethod::GetMaxNumPaces()
{
    return mMaxNumPaces;
}

void AbstractActionPotentialMethod::Reset() { mRunYet = false; }

void AbstractActionPotentialMethod::SetLackOfOneToOneCorrespondenceIsError(
    bool errorOn)
{
    mNoOneToOneCorrespondenceIsError = errorOn;
}

void AbstractActionPotentialMethod::SetAlternansIsError(bool errorOn)
{
    mAlternansIsError = errorOn;
}

bool AbstractActionPotentialMethod::DidErrorOccur(void)
{
    if (!mRunYet)
    {
        EXCEPTION("Simulation has not been run yet.");
    }

    return !mSuccessful;
}

std::string AbstractActionPotentialMethod::GetErrorMessage()
{
    if (!mRunYet)
    {
        EXCEPTION("Simulation has not been run yet.");
    }

    if (mSuccessful)
    {
        EXCEPTION("No error occurred.");
    }
    return mErrorMessage;
}

unsigned AbstractActionPotentialMethod::GetErrorCode()
{
    if (!mRunYet)
    {
        EXCEPTION("Simulation has not been run yet.");
    }

    if (mSuccessful)
    {
        return 0u;
    }
    return mErrorCode;
}

void AbstractActionPotentialMethod::WriteMessageToFile(
    const std::string& rMessage)
{
    WARNING(rMessage); // Send the message to std::cout as well as any message
    // file in subclasses.
}

void AbstractActionPotentialMethod::SuppressOutput(bool suppress)
{
    mSuppressOutput = suppress;
}

OdeSolution AbstractActionPotentialMethod::SteadyStatePacingExperiment(
    boost::shared_ptr<AbstractCvodeCell> pModel, double& rApd90, double& rApd50,
    double& rUpstroke, double& rPeak, double& rPeakTime, double& rCaMax,
    double& rCaMin, const double printingTimeStep, const double conc)
{
    mRunYet = true;
    mRepeat = false; //    These two resets are for counting whether we should repeat the last
    mRepeatNumber = 0u; // AP simulation to detect alternans or abnormal repolarization.

    /**
     * STEADY STATE PACING EXPERIMENT
     *
     * Stimulate many times to establish the steady state response.
     * Do most of the calculations in one go to avoid overhead setting up Cvode...
     */
    bool skip_steady_state_calculation = false;

    // Get the stimulus parameters and make sure the maximum timestep is the
    // minimum of
    // * stimulus duration
    // or
    // * printing time step.

    if (!boost::dynamic_pointer_cast<RegularStimulus>(
            pModel->GetStimulusFunction()))
    {
        EXCEPTION(
            "AbstractActionPotentialMethod only works with cells that have a "
            "RegularStimulus set.");
    }

    boost::shared_ptr<RegularStimulus> p_default_stimulus = boost::static_pointer_cast<RegularStimulus>(
        pModel->GetStimulusFunction());
    unsigned voltage_index = pModel->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");

    double s1_period = p_default_stimulus->GetPeriod();
    double s_duration = p_default_stimulus->GetDuration();
    double maximum_time_step; // ms
    if (printingTimeStep > s_duration)
    {
        maximum_time_step = s_duration;
    }
    else
    {
        maximum_time_step = printingTimeStep;
    }

    const unsigned num_paces_to_analyze = 2u;

    // In this block we just do one/two complete action potentials
    // to see if the model 'fires off' an action potential properly.
    // If it doesn't we don't bother pursuing this concentration and we skip
    // running to steady state.
    {
        pModel->SetMaxSteps(100000);
        OdeSolution solution = pModel->Solve(0, num_paces_to_analyze * s1_period, maximum_time_step,
                                             printingTimeStep); // Maximum timestep here is usually
        // the printing time step
        std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);

        // Get voltage properties using an action potential threshold
        // If this is the control drug case set a sensible threshold for APs,
        // otherwise use pre-existing, or default (-50mV as set in constructor).
        if (fabs(conc) < 1e-10 && !mActionPotentialThresholdSetManually)
        {
            double max_v = -DBL_MAX;
            double min_v = DBL_MAX;
            for (unsigned i = 0; i < voltages.size(); i++)
            {
                if (voltages[i] > max_v)
                {
                    max_v = voltages[i];
                }
                if (voltages[i] < min_v)
                {
                    min_v = voltages[i];
                }
            }
            mActionPotentialThreshold = min_v + (max_v - min_v) * 0.1; // Say threshold for an AP is 10% of way up upstroke.
        }

        CellProperties voltage_properties(voltages, solution.rGetTimes(),
                                          mActionPotentialThreshold);
        double apd90;
        try
        {
            apd90 = voltage_properties.GetLastActionPotentialDuration(90);
            if (!mSuppressOutput)
                std::cout << "First pace APD90 = " << apd90 << "\n"; // << std::flush;
        }
        catch (Exception& e)
        {
            if (e.GetShortMessage() == "AP did not occur, never exceeded threshold voltage." || e.GetShortMessage() == "No full action potential was recorded")
            {
                skip_steady_state_calculation = true;
            }
            else
            {
                throw e;
            }
        }
    }

    // The mMaxNumPaces is based on the total we want to do, so the steady state
    // runner
    // should go for fewer, as some are analysed before this, and some after.
    const unsigned num_paces_analysed_elsewhere = 2u * num_paces_to_analyze;

    if (!skip_steady_state_calculation && mMaxNumPaces > num_paces_analysed_elsewhere)
    {
        // This method tries to detect a steady state if we are happy we're
        // producing APs...
        SteadyStateRunner steady_runner(pModel, true); // The true denotes a 'two
        // pace scan', which will
        // catch alternans too.
        if (mSuppressOutput)
        {
            steady_runner.SuppressOutput();
        }
        if (mMaxNumPaces != UNSIGNED_UNSET)
        {
            steady_runner.SetMaxNumPaces(mMaxNumPaces - num_paces_analysed_elsewhere);
        }
        steady_runner.RunToSteadyState();
    }

    OdeSolution solution = PerformAnalysisOfTwoPaces(
        pModel, rApd90, rApd50, rUpstroke, rPeak, rPeakTime, rCaMax, rCaMin,
        s1_period, maximum_time_step, printingTimeStep, conc);

    if (mRepeat)
    {
        // std::cout << "Repeating simulation to order alternans APs consistently...\n";
        // If we might benefit from pushing forward one period and re-analysing...
        PushModelForwardOneS1Interval(pModel, s1_period, maximum_time_step);

        solution = PerformAnalysisOfTwoPaces(
            pModel, rApd90, rApd50, rUpstroke, rPeak, rPeakTime, rCaMax, rCaMin,
            s1_period, maximum_time_step, printingTimeStep, conc);
    }

    return solution;
}

void AbstractActionPotentialMethod::PushModelForwardOneS1Interval(
    boost::shared_ptr<AbstractCvodeCell> pModel, double pacingCycleLength,
    double maxTimeStep)
{
    pModel->Solve(0, pacingCycleLength, maxTimeStep);
}

OdeSolution AbstractActionPotentialMethod::PerformAnalysisOfTwoPaces(
    boost::shared_ptr<AbstractCvodeCell> pModel, double& rApd90, double& rApd50,
    double& rUpstroke, double& rPeak, double& rPeakTime, double& rCaMax,
    double& rCaMin, const double s1_period, const double maximumTimeStep,
    const double printingTimeStep, const double conc)
{
    // We usually analyse two paces to look for alternans.
    // Have observed three period behaviour or more, so this is a hardcoded option
    // for now.
    const unsigned num_paces_to_analyze = 2u;
    mRepeat = false;
    const double alternans_threshold = 1; // ms in APD90 - hardcoded,
    // could make an option in future. But if it is any smaller alternans 'comes
    // and goes' as you
    // move through parameter space. Here it appears to pick up only the serious
    // (after a bifurcation)
    // kind of alternans.

    pModel->SetMaxSteps(num_paces_to_analyze * 100000); // Internal ODE solver steps, not paces!
    OdeSolution solution = pModel->Solve(0, num_paces_to_analyze * s1_period, maximumTimeStep,
                                         printingTimeStep); // Get plenty of detail on these two
    // paces for analysis.

    // Get voltage properties
    const unsigned voltage_index = pModel->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
    std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
    CellProperties voltage_properties(voltages, solution.rGetTimes(),
                                      mActionPotentialThreshold);

    std::vector<double> apd90s;
    std::vector<double> peak_voltages;
    // See if we can get back some action potential duration(s).
    try
    {
        apd90s = voltage_properties.GetAllActionPotentialDurations(90);
        if (!mSuppressOutput)
        {
            std::cout << "Last " << apd90s.size()
                      << " AP duration(s) = "; //<< std::flush;
            for (unsigned i = 0; i < apd90s.size(); i++)
            {
                std::cout << apd90s[i] << ",\t";
            }
            std::cout << std::endl; //<< std::flush;
        }

        if (apd90s.size() >= 2u && fabs(apd90s[0] - apd90s[1]) > alternans_threshold)
        {
            // We suspect alternans, and analyse the first of the two APs
            rApd90 = voltage_properties.GetAllActionPotentialDurations(90)[0];
            rApd50 = voltage_properties.GetAllActionPotentialDurations(50)[0];
            rUpstroke = voltage_properties.GetMaxUpstrokeVelocities()[0];
            peak_voltages = voltage_properties.GetPeakPotentials();
            rPeak = peak_voltages[0];
            rPeakTime = voltage_properties.GetTimesAtPeakPotentials()[0];
        }
        else
        {
            // Return the last as it is more likely to be the steady state one.
            rApd90 = voltage_properties.GetLastActionPotentialDuration(90);
            rApd50 = voltage_properties.GetLastActionPotentialDuration(50);
            rUpstroke = voltage_properties.GetLastCompleteMaxUpstrokeVelocity();
            rPeak = voltage_properties.GetLastCompletePeakPotential();
            rPeakTime = voltage_properties.GetTimeAtLastCompletePeakPotential();
        }
        // It makes sense to return the peak voltage time relative to start of
        // stimulus application.
        boost::shared_ptr<RegularStimulus> p_reg_stim = boost::static_pointer_cast<RegularStimulus>(
            pModel->GetStimulusFunction());
        rPeakTime = std::fmod(rPeakTime - p_reg_stim->GetStartTime(), s1_period);

        if (pModel->HasAnyVariable("cytosolic_calcium_concentration"))
        {
            if (pModel->HasDerivedQuantity("cytosolic_calcium_concentration"))
            {
                solution.CalculateDerivedQuantitiesAndParameters(pModel.get());
            }
            std::vector<double> calcium = solution.GetAnyVariable("cytosolic_calcium_concentration");
            rCaMax = *(std::max_element(calcium.begin(), calcium.end()));
            rCaMin = *(std::min_element(calcium.begin(), calcium.end()));
        }
        else
        {
            WARN_ONCE_ONLY(pModel->GetSystemName()
                           << " does not have 'cytosolic_calcium_concentration' "
                              "annotated, please tag it if it is present.");
        }
        mSuccessful = true;
    }
    catch (Exception& e)
    {
        // We didn't get back any action potentials,
        // which should be because of one of the following errors (if not throw the
        // exception).
        if (e.GetShortMessage() != "AP did not occur, never exceeded threshold voltage."
            && e.GetShortMessage() != "No full action potential was recorded"
            && e.GetShortMessage() != "No MaxUpstrokeVelocity matching a full action potential was recorded.")
        {
            throw e;
        }

        std::stringstream message;
        if (conc != DOUBLE_UNSET)
        {
            message << "At a concentration of " << conc << "uM: ";
        }
        message << "no action potentials were recorded, cell did not ";

        // Work out whether most of the time was spent above or below threshold
        double mean_voltage = 0.0;
        for (unsigned i = 0; i < voltages.size(); i++)
        {
            mean_voltage += voltages[i];
        }
        mean_voltage /= ((double)(voltages.size()));

        if (mean_voltage > mActionPotentialThreshold)
        {
            mErrorCode = 2u;
            mErrorMessage = "NoActionPotential_2";
            message << "repolarise.";
        }
        else
        {
            mErrorCode = 1u;
            mErrorMessage = "NoActionPotential_1";
            message << "depolarise.";
        }
        // std::cout << message.str() << std::endl << std::flush;
        if (!mSuppressOutput)
            std::cout << message.str() << std::endl
                      << std::flush;
        WriteMessageToFile(message.str());
        mSuccessful = false;
        mPeriodTwoBehaviour = true;
        // Redo the analysis in case being 'in sync' on period 2 orbit allows us
        // to get an APD.
        mRepeat = true;
    }

    if (mSuccessful)
    {
        assert(apd90s.size() > 0u);

        // Deal with the case when we are in alternans
        if (apd90s.size() >= 2u && fabs(apd90s[0] - apd90s[1]) > alternans_threshold)
        {
            mPeriodTwoBehaviour = true;
            if (apd90s[1] > apd90s[0] && mRepeatNumber == 0u)
            {
                // Redo so that we always plot the longest AP first.
                mRepeat = true;
            }
            else
            { // If we're going to repeat, we don't want this message twice.
                if (mAlternansIsError)
                {
                    // These conditions check (if we are looking for repolarisation caused
                    // alternans by setting
                    // a default APD90):
                    // 1. That the alternans APDs are longer than default
                    // 2. That this isn't caused by a slow depolarisation which tends to
                    // give varied Peak Vm
                    // (see test case 9 in O'Hara model TestTroublesomeApEvaluations)
                    if (mDefaultParametersApd90 != DOUBLE_UNSET && (apd90s[0] > mDefaultParametersApd90 && apd90s[1] > mDefaultParametersApd90) && (fabs(peak_voltages[0] - peak_voltages[1]) < 10 /*mV*/)) // magic number!
                    {
                        // We have alternans tending to long/no repolarisation
                        mErrorCode = 6u;
                        mErrorMessage = "NoActionPotential_6";
                    }
                    else
                    {
                        mErrorCode = 4u;
                        mErrorMessage = "NoActionPotential_4";
                    }
                    mSuccessful = false;
                }

                std::stringstream message;
                if (conc != DOUBLE_UNSET)
                {
                    message << "At a concentration of " << conc << "uM: ";
                }
                message << "possible alternans detected, APD90s = " << apd90s[0] << ", "
                        << apd90s[1] << " ms";
                std::string message_string = message.str();
                WriteMessageToFile(message_string);
            }

            // If we want to create an error code for 'struggling to depolarise'
            // This will overwrite alternans errors (it's more useful).
            if (mDefaultParametersTimeOfVMax != DOUBLE_UNSET && rPeakTime > mDefaultParametersTimeOfVMax + 80) // ms MAGIC NUMBER!
            {
                mErrorCode = 7u;
                mErrorMessage = "NoActionPotential_7";
                mSuccessful = false;
            }
        }

        // Deal with the case when there was only one AP but two stimuli
        if (apd90s.size() < num_paces_to_analyze && mSuccessful)
        {
            std::stringstream message;
            if (conc != DOUBLE_UNSET)
            {
                message << "At a concentration of " << conc << "uM: ";
            }
            message << "only one action potential was recorded (" << apd90s[0]
                    << " ms) for two paces of " << s1_period << "ms.";
            std::string message_string = message.str();
            WriteMessageToFile(message_string);

            if (mNoOneToOneCorrespondenceIsError)
            {
                // The second condition here is to check if the second AP is
                // non-repolarising...
                // absence of this caused a bit of a bug where a period 3
                // non-repolarising looked like non
                // depolarising...
                if (apd90s[0] > s1_period || voltages.back() >= mActionPotentialThreshold)
                {
                    mErrorCode = 3u;
                    mErrorMessage = "NoActionPotential_3";
                }
                else
                {
                    mErrorCode = 5u;
                    mErrorMessage = "NoActionPotential_5";
                }
                mSuccessful = false;
            }
            mPeriodTwoBehaviour = true;
        }
    }

    mRepeatNumber++;

    return solution;
}

void AbstractActionPotentialMethod::
    SetVoltageThresholdForRecordingAsActionPotential(double threshold)
{
    mActionPotentialThreshold = threshold;
    mActionPotentialThresholdSetManually = true;
}

void AbstractActionPotentialMethod::SetControlActionPotentialDuration90(
    double apd90)
{
    mDefaultParametersApd90 = apd90;
}

void AbstractActionPotentialMethod::SetControlTimeOfPeakVoltage(
    double timeVMax)
{
    mDefaultParametersTimeOfVMax = timeVMax;
}
