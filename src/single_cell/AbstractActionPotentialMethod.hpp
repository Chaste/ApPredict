/*

Copyright (c) 2005-2021, University of Oxford.
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

#ifndef ABSTRACTACTIONPOTENTIALMETHOD_HPP_
#define ABSTRACTACTIONPOTENTIALMETHOD_HPP_

#include <boost/shared_ptr.hpp>
#include <vector>

#include "AbstractCvodeCell.hpp"
#include "CellProperties.hpp"
#include "CommandLineArguments.hpp"
#include "Exception.hpp"
#include "OdeSolution.hpp"
#include "RegularStimulus.hpp"
#include "SteadyStateRunner.hpp"

/**
 * This is a class that provides methods for running cell models to steady state
 * and evaluating
 * the resulting action potential markers like APD90.
 */
class AbstractActionPotentialMethod
{
private:
    /** Whether we have run a simulation */
    bool mRunYet;

    /**
     * The limit on the maximum number of paces to take
     *
     * Value is UNSIGNED_UNSET if this has not been set.
     */
    unsigned mMaxNumPaces;

    /** Any error messages that resulted */
    std::string mErrorMessage;

    /** The error code (same as end of mErrorMessage) or 0 if no error occurred */
    unsigned mErrorCode;

    /** Whether to say that no 1:1 correspondence between Stimuli and APs is an
     * error, default is no, just warning */
    bool mNoOneToOneCorrespondenceIsError;

    /** Whether to say that alternans behaviour is an error code. Default is no,
     * just a warning. */
    bool mAlternansIsError;

    /** The threshold for an action potential (defaults to -50mV), altered if
     * model is run at control conc of drug */
    double mActionPotentialThreshold;

    /** Whether the AP threshold was set manually (defaults to false) */
    bool mActionPotentialThresholdSetManually;

    /** The default APD90, can optionally be set to get more accurate error
     * messages about whether
     * errors are to do with depolarisation or repolarisation.
     */
    double mDefaultParametersApd90;

    /** The default dVdt_Max, can optionally be set to get more accurate error
     * messages and warnings. */
    double mDefaultParametersTimeOfVMax;

    /** Whether to reanalyse one period further on if in alternans */
    bool mRepeat;

    /** We attempt to reanalyse one period further on if in alternans */
    unsigned mRepeatNumber;

    /**
     * A helper method, only available to this class.
     *
     * @param pModel  A boost shared pointer to a cardiac cell model, with a
     * RegularStimulus
     * @param rApd90  double to populate with the postprocessed APD90 (ms)
     * @param rApd50  double to populate with the postprocessed APD50 (ms)
     * @param rUpstroke  double to populate with the postprocessed maximum
     * upstroke velocity (mV/ms)
     * @param rPeak  double to populate with the postprocessed peak voltage (mV)
     * @param rPeakTime  double to populate with the time of postprocessed peak voltage (mV)
     * @param rCaMax double to populate with the maximum cytosolic calcium
     * concentration
     * @param rCaMin double to populate with the minimum cytosolic calcium
     * concentration
     * @param maximumTimeStep  The maximum CVODE time step to use (ms).
     * @param printingTimeStep  the printing time step to use (defaults to 1ms).
     * @param conc  [optional] concentration argument (only used for more helpful
     * warning messages).
     *
     * @return the solution of the ODE.
     */
    OdeSolution PerformAnalysisOfTwoPaces(
        boost::shared_ptr<AbstractCvodeCell> pModel, double& rApd90,
        double& rApd50, double& rUpstroke, double& rPeak, double& rPeakTime, double& rCaMax,
        double& rCaMin, const double s1_period, const double maximumTimeStep,
        const double printingTimeStep, const double conc);

    /**
     * A little method to 'push' cell model forward one S1 period, to get it 'in
     * sync'
     * if it is doing different alternans, or capturing the end of an AP in 2:1
     * pacing:AP.
     *
     * @param pModel  A boost shared pointer to a cardiac cell model, with a
     * RegularStimulus
     * @param pacingCycleLength  The period of the regular stimulus in this model
     * @param maxTimeStep  The maximum CVODE time step to use in solving
     */
    void PushModelForwardOneS1Interval(
        boost::shared_ptr<AbstractCvodeCell> pModel, double pacingCycleLength,
        double maxTimeStep);

protected:
    /** Whether to suppress output to std::cout */
    bool mSuppressOutput;

    /** The frequency in Hz at which to perform this run */
    double mHertz;

    /** Whether the run was successful */
    bool mSuccessful;

    /** Whether anything happened on a period-2 basis rather than period 1, e.g.
     * alternans, long APs etc. */
    bool mPeriodTwoBehaviour;

    /**
     * This method conducts a steady-state single cell simulation on a model
     * (using its internal stimulus,
     * which must be a RegularStimulus).
     *
     * It populates a number of steady state action potential property arguments
     * and puts warnings to screen, or file
     * if specified by subclasses.
     *
     * The mErrorMessage takes one of the following strings:
     *  * NoActionPotential_1 - cell did not depolarise enough to form action potential
     *
     *  * NoActionPotential_2 - cell did not repolarise enough to form action potential
     *
     *  * NoActionPotential_3 - action potential did not repolarise in time for the next stimulus (no 1:1 correspondence) only an error if
     * #SetLackOfOneToOneCorrespondenceIsError() has been called, otherwise it's a warning.
     *
     *  * NoActionPotential_4 - alternans or period two behaviour is occurring. This
     * is only an error if SetAlternansIsError() has been called, otherwise it is a warning. If
     * SetControlActionPotentialDuration90 has been set, then this code indicates short alternans and
     * NoActionPotential_6 indicates long alternans.
     *
     * * NoActionPotential_5 - alternans or period two behaviour is occurring, but
     * one AP is sub-threshold. This is only an error if #SetLackOfOneToOneCorrespondenceIsError() has been
     * called,otherwise it is a warning.
     *
     * * NoActionPotential_6 - long alternans or period two behaviour is occurring.
     * This is only an error if #SetLackOfOneToOneCorrespondenceIsError() has been
     * called, otherwise it is a warning. You only get this error code if
     * SetControlActionPotentialDuration90() has been
     * called, otherwise all alternans will be NoActionPotential_4.
     *
     * * NoActionPotential_7 -  if depolarisation abnormalities are occurring
     * You will only get this code if SetControlTimeOfPeakVoltage()
     * has been set.
     *
     * Optionally you can state a concentration (for more helpful warning
     * messages) and printing time step.
     *
     * @param pModel  A boost shared pointer to a cardiac cell model, with a RegularStimulus.
     * @param rApd90  double to populate with the postprocessed APD90 (ms)
     * @param rApd50  double to populate with the postprocessed APD50 (ms)
     * @param rUpstroke  double to populate with the postprocessed maximum upstroke velocity (mV/ms)
     * @param rPeak  double to populate with the postprocessed peak voltage (mV)
     * @param rPeak  double to populate with the time of postprocessed peak voltage (mV)
     * @param rCaMax  double to populate with the maximum of calcium transient (mM)
     * @param rCaMin  double to populate with the minimum of calcium transient (mM)
     * @param printingTimeStep  the printing time step to use (defaults to 1ms).
     * @param conc  [optional] concentration argument (only used for more helpful warning messages).
     *
     * @return the solution of the ODE.
     */
    OdeSolution SteadyStatePacingExperiment(
        boost::shared_ptr<AbstractCvodeCell> pModel, double& rApd90,
        double& rApd50, double& rUpstroke, double& rPeak, double& rPeakTime, double& rCaMax,
        double& rCaMin,
        const double printingTimeStep = 1, // ms
        const double conc = DOUBLE_UNSET);

    /**
      * This method just passes any message into the WARNINGS macro.
      *
      * It can be called from subclass methods to write a log message to the
      * messages.txt file
      * that should be displayed alongside the results.
      * (for example a warning that the cell failed to de/re-polarise at a certain
      * concentration.)
      * @param rMessage  The message to write to file
      */
    virtual void WriteMessageToFile(const std::string& rMessage);

public:
    /**
    * The main constructor, sets a few defaults and reads some command line
    * options.
    */
    AbstractActionPotentialMethod();

    /**
    * Destructor
    */
    virtual ~AbstractActionPotentialMethod(){};

    /**
    * Set the maximum number of paces to run for.
    * Use of this method will override any command line modifications
    * set in the constructor.
    *
    * @param numPaces  The maximum number of paces to run in an attempt to get to
    * steady state.
    */
    void SetMaxNumPaces(unsigned numPaces);

    /**
    * @return the maximum number of paces to run to steady state.
    */
    unsigned GetMaxNumPaces();

    /**
    * @return  whether an error occurred in the action potential marker
    * evaluation (not simulation itself).
    */
    bool DidErrorOccur(void);

    /**
     * @return The error message from attempting to evaluate markers, such as
     * APD90.
     *
     * This is "NoActionPotential_1" if the cell failed to depolarise (no AP
     * occurs).
     *
     * Or "NoActionPotential_2" if the cell failed to repolarise (AP never returns
     * to resting potential).
     *
     * Or "NoActionPotential_3" action potential did not repolarise in time for
     * the next stimulus (no 1:1 correspondence).
     *         This only occurs if #SetLackOfOneToOneCorrespondenceIsError() has
     * been called.
     *
     * Or "NoActionPotential_4" if alternans is occurring (different APD on
     * subsequent beats).
     * This only occurs if #SetAlternansIsError() has been called,
     * otherwise a warning is printed to screen.
     * If SetControlActionPotentialDuration90() has been set, then this code
     * indicates short
     * alternans and NoActionPotential_6 indicates long alternans. Otherwise this
     * code covers both
     *
     * Or "NoActionPotential_5" if alternans is occurring but second AP is
     * sub-threshold on
     * the next stimulus (no 1:1 correspondence).
     * This only occurs if #SetLackOfOneToOneCorrespondenceIsError() has been
     * called.
     *
     * Or "NoActionPotential_6" - long alternans or period two behaviour is
     * occurring.
     * This is only an error if #SetLackOfOneToOneCorrespondenceIsError() has been
     * called,
     * otherwise it is a warning. You will only get this code if
     * SetControlActionPotentialDuration90()
     * has been set.
     *
     * Or "NoActionPotential_7" if depolarisation abnormalities are occurring
     * You will only get this code if SetControlTimeOfPeakVoltage()
     * has been set.
     */
    std::string GetErrorMessage();

    /**
     * @return The error code from attempting to evaluate markers, such as APD90.
     *
     * Returns "0" if no error occurred.
     *
     * This is "1" if the cell failed to depolarise (no AP occurs).
     *
     * Or "2" if the cell failed to repolarise (AP never returns to resting
     * potential).
     *
     * Or "3" action potential did not repolarise in time for the next stimulus
     * (no 1:1 correspondence).
     *         This only occurs if #SetLackOfOneToOneCorrespondenceIsError() has
     * been called.
     *
     * Or "4" if alternans is occurring (different APD on subsequent beats).
     * This only occurs if #SetAlternansIsError() has been called,
     * otherwise a warning is printed to screen. If
     * SetControlActionPotentialDuration90()
     * has been set you get this error code for shorter than average APDs, if not
     * you
     * get this error code for short and long alternans.
     *
     * Or "5" if alternans is occurring but second AP is sub-threshold on
     * the next stimulus (no 1:1 correspondence).
     * This only occurs if #SetLackOfOneToOneCorrespondenceIsError() has been
     * called.
     *
     * Or "6" if alternans is occurring (with longer than default APD90s)
     * This only occurs if #SetLackOfOneToOneCorrespondenceIsError() has been
     * called.
     * You will only get this code if SetControlActionPotentialDuration90()
     * has been set.
     *
     * Or "7" if depolarisation abnormalities are occurring
     * You will only get this code if SetControlTimeOfPeakVoltage()
     * has been set.
     */
    unsigned GetErrorCode();

    /**
     * Set the voltage which is considered to be a triggered action potential.
     *
     * If this is not set it defaults to -50mV.
     *
     * @param threshold
     */
    void SetVoltageThresholdForRecordingAsActionPotential(double threshold);

    /**
     * Set the control APD90. This is used in parameter sweeps to decide whether an
     * error message should be more like "no repolarisation" or "no
     * depolarisation"!
     *
     * If this is not set, it is not used.
     *
     * @param apd90  Default parameters APD90.
     */
    void SetControlActionPotentialDuration90(double apd90);

    /**
     * Set the control dV/dt max. This is used in parameter sweeps to decide whether an
     * error message should be more like "no repolarisation" or "no
     * depolarisation"!
     *
     * If this is not set, it is not used.
     *
     * @param dvDtMax  max(dV/dt) with default parameters.
     */
    void SetControlTimeOfPeakVoltage(double dvDtMax);

    /**
     * Tell the simulator whether to suppress output to std::cout.
     *
     * @param suppress  Whether to suppress output (defaults to true).
     */
    void SuppressOutput(bool suppress = true);

    /**
     * Tell this class to treat a lack of 1:1 correspondence between stimuli and
     * action potentials
     * as an error, rather than just a warning.
     */
    void SetLackOfOneToOneCorrespondenceIsError(bool errorOn = true);

    /**
     * Tell this class to treat alternans / period two behaviour as an error,
     * rather than just a warning.
     */
    void SetAlternansIsError(bool errorOn = true);

    /**
     * Reset the action potential evaluator for a new run.
     */
    void Reset();
};

#endif // ABSTRACTACTIONPOTENTIALMETHOD_HPP_
