/*

Copyright (c) 2005-2013, University of Oxford.
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

#include "AbstractActionPotentialMethod.hpp"

AbstractActionPotentialMethod::AbstractActionPotentialMethod()
     : mMaxNumPacesSet(false),
       mRunYet(false),
       mMaxNumPaces(UNSIGNED_UNSET),
       mErrorMessage(""),
       mNoOneToOneCorrespondenceIsError(false),
       mActionPotentialThreshold(-50),
       mSuppressOutput(false),
       mHertz(1.0), // Hz
       mSuccessful(false)
{
    CommandLineArguments* p_args = CommandLineArguments::Instance();

    if (p_args->OptionExists("--pacing-freq"))
    {
        mHertz = p_args->GetDoubleCorrespondingToOption("--pacing-freq");
    }
    //if (!mSuppressOutput) std::cout << "* Pacing Frequency = " << mHertz << " Hz\n";

    if (p_args->OptionExists("--pacing-max-time"))
    {
        // Max number of paces is equal to maximum pacing time (converted from minutes to seconds) x paces per second.
        SetMaxNumPaces((unsigned)(60.0*mHertz*p_args->GetDoubleCorrespondingToOption("--pacing-max-time")));
    }
}

void AbstractActionPotentialMethod::SetMaxNumPaces(unsigned numPaces)
{
    if (numPaces != UNSIGNED_UNSET)
    {
        mMaxNumPaces = numPaces;
        mMaxNumPacesSet = true;
    }
}

void AbstractActionPotentialMethod::Reset()
{
    mRunYet = false;
}

void AbstractActionPotentialMethod::SetLackOfOneToOneCorrespondenceIsError(bool errorOn)
{
    mNoOneToOneCorrespondenceIsError = errorOn;
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

void AbstractActionPotentialMethod::WriteMessageToFile(const std::string& rMessage)
{
    WARNING(rMessage); // Send the message to std::cout as well as any message file in subclasses.
}

void AbstractActionPotentialMethod::SuppressOutput(bool suppress)
{
    mSuppressOutput = suppress;
}


OdeSolution AbstractActionPotentialMethod::SteadyStatePacingExperiment(
                                     boost::shared_ptr<AbstractCvodeCell> pModel,
                                     double& rApd90,
                                     double& rApd50,
                                     double& rUpstroke,
                                     double& rPeak,
                                     const double printingTimeStep,
                                     const double conc)
{
    mRunYet = true;

    /**
     * STEADY STATE PACING EXPERIMENT
     *
     * Stimulate many times to establish the steady state response.
     * Do most of the calculations in one go to avoid overhead setting up Cvode...
     */
    bool skip_steady_state_calculation = false;

    // Get the stimulus parameters and make sure the maximum timestep is the minimum of
    // * stimulus duration
    // or
    // * printing time step.

    if (!boost::dynamic_pointer_cast<RegularStimulus>(pModel->GetStimulusFunction()))
    {
        EXCEPTION("AbstractActionPotentialMethod only works with cells that have a RegularStimulus set.");
    }

    boost::shared_ptr<RegularStimulus> p_default_stimulus =
                         boost::static_pointer_cast<RegularStimulus>(pModel->GetStimulusFunction());
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

    // In this block we just do one complete action potential
    // to see if the model 'fires off' an action potential properly.
    // If it doesn't we don't bother pursuing this concentration and we skip running to steady state.
    {
        pModel->SetMaxSteps(100000);
        OdeSolution solution = pModel->Solve(0,s1_period, maximum_time_step, printingTimeStep); // Maximum timestep here is usually the printing time step
        std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);

        // Get voltage properties using an action potential threshold
        // If this is the control drug case set a sensible threshold for APs,
        // otherwise use pre-existing, or default (-50mV as set in constructor).
        if (fabs(conc)<1e-10)
        {
            double max_v = -DBL_MAX;
            double min_v = DBL_MAX;
            for (unsigned i=0 ; i<voltages.size(); i++)
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
            mActionPotentialThreshold = min_v + (max_v - min_v)*0.1; // Say threshold for an AP is 10% of way up upstroke.
        }

        CellProperties voltage_properties(voltages, solution.rGetTimes(), mActionPotentialThreshold);
        double apd90;
        try
        {
            apd90 = voltage_properties.GetLastActionPotentialDuration(90);
            if (!mSuppressOutput) std::cout << "First pace APD90 = " << apd90 << "\n";// << std::flush;
        }
        catch (Exception& e)
        {
            if (e.GetShortMessage()=="AP did not occur, never exceeded threshold voltage." ||
                e.GetShortMessage()=="No full action potential was recorded")
            {
                skip_steady_state_calculation = true;
                //EXCEPTION("No action potential on first pace - model/parameter choice is dodgy\n");
            }
            else
            {
                throw e;
            }
        }
    }

    // This method tries to detect a steady state if we are happy we're producing APs...
    if (!skip_steady_state_calculation)
    {
        SteadyStateRunner steady_runner(pModel, true); // The true denotes a 'two pace scan' which will catch alternans too.
        if (mMaxNumPacesSet)
        {
            steady_runner.SetMaxNumPaces(mMaxNumPaces);
        }
        steady_runner.RunToSteadyState();
    }

    // Set the max time step to be the printing time step (1ms).
    // Do one final pace that we will analyse...
    const unsigned num_paces_to_analyze = 2u;
    pModel->SetMaxSteps(num_paces_to_analyze*100000);
    OdeSolution solution = pModel->Solve(0,num_paces_to_analyze*s1_period, maximum_time_step, printingTimeStep); // Get plenty of detail on these two paces for analysis.

    // Get voltage properties
    double threshold = -50; // Threshold for AP is above -50ish.
    std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
    CellProperties voltage_properties(voltages, solution.rGetTimes(), threshold);

    std::vector<double> apd90s;
    try
    {
        apd90s = voltage_properties.GetAllActionPotentialDurations(90);
        if (!mSuppressOutput)
        {
            std::cout << "Last " << apd90s.size() << " AP duration(s) = ";//<< std::flush;
            for (unsigned i=0; i<apd90s.size(); i++)
            {
                std::cout << apd90s[i] << ",\t";
            }
            std::cout << std::endl ;//<< std::flush;
        }
        //apd50_0 = voltage_properties.GetAllActionPotentialDurations(50)[0];
        //upstroke_0 = voltage_properties.GetAllCompleteMaxUpstrokeVelocities()[0];
        //peak_0 = voltage_properties.GetAllCompletePeakPotentials()[0];
        rApd90 = voltage_properties.GetLastActionPotentialDuration(90);
        rApd50 = voltage_properties.GetLastActionPotentialDuration(50);
        rUpstroke = voltage_properties.GetLastCompleteMaxUpstrokeVelocity();
        rPeak = voltage_properties.GetLastCompletePeakPotential();
        mSuccessful = true;
    }
    catch (Exception& e)
    {
        if (e.GetShortMessage()=="AP did not occur, never exceeded threshold voltage." ||
            e.GetShortMessage()=="No full action potential was recorded" ||
            e.GetShortMessage()=="No MaxUpstrokeVelocity matching a full action potential was recorded.")
        {
            std::stringstream message;
            if (conc!=DOUBLE_UNSET)
            {
                message << "At a concentration of " << conc << "uM: ";
            }
            message << "no action potentials were recorded, cell did not ";

            // Work out whether most of the time was spent above or below threshold
            double mean_voltage = 0.0;
            for (unsigned i=0; i<voltages.size(); i++)
            {
                mean_voltage += voltages[i];
            }
            mean_voltage /= ((double)(voltages.size()));

            if (mean_voltage > threshold)
            {
                mErrorMessage = "NoActionPotential_2";
                message << "repolarise.";
            }
            else
            {
                mErrorMessage = "NoActionPotential_1";
                message << "depolarise.";
            }
            //std::cout << message.str() << std::endl << std::flush;
            if (!mSuppressOutput) std::cout << message.str() << std::endl << std::flush;
            WriteMessageToFile(message.str());
            mSuccessful = false;
        }
        else
        {
            throw e;
        }
    }

    if (mSuccessful)
    {
        assert(apd90s.size()>0u);
        if (apd90s.size() >= 2u && fabs(apd90s[0]-apd90s[1])>0.5)
        {
            std::stringstream message;
            if (conc!=DOUBLE_UNSET)
            {
                message << "At a concentration of " << conc << "uM:";
            }
            message << "possible alternans detected, APD90s = " << apd90s[0] << ", " << apd90s[1] << " ms";
            std::string message_string = message.str();
            WriteMessageToFile(message_string);
        }

        if (apd90s.size() < num_paces_to_analyze)
        {
            std::stringstream message;
            if (conc!=DOUBLE_UNSET)
            {
                message << "At a concentration of " << conc << "uM:";
            }
            message << "only one action potential was recorded (" << apd90s[0] << " ms) for two paces of " << s1_period << "ms.";
            std::string message_string = message.str();
            WriteMessageToFile(message_string);

            if (mNoOneToOneCorrespondenceIsError)
            {
                mErrorMessage = "NoActionPotential_3";
                mSuccessful = false;
            }
        }
    }

    return solution;
}

