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

#include <fstream>
#include "ApPredictMethods.hpp"

#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "CommandLineArguments.hpp"
#include "FileFinder.hpp"

#include "AbstractDataStructure.hpp"
#include "DoseCalculator.hpp"
#include "ActionPotentialDownsampler.hpp"

#include "ZeroStimulus.hpp"
#include "RegularStimulus.hpp"

#include "SetupModel.hpp"
#include "SteadyStateRunner.hpp"
#include "ProgressReporter.hpp"


std::string ApPredictMethods::PrintArguments()
{
	std::string message = "ApPredict::Please provide these inputs:\n" + SetupModel::PrintArguments();
	message += PrintCommonArguments();
	return message;
}

std::string ApPredictMethods::PrintCommonArguments()
{
    std::string message = "*\n"
            "* Specifying pacing frequency:\n"
            "* --pacing-freq       Pacing frequency (Hz) (optional - defaults to 1Hz)\n"
            "* --pacing-max-time   Maximum time for which to pace the cell model in MINUTES (optional)\n"
            "*\n"
            "* Specifying drug/ion-channel dose-response properties for each channel:\n"
            "* Channels are named:\n"
            "* * herg (IKr current - hERG),\n"
            "* * na (fast sodium current - NaV1.5),\n"
            "* * cal (L-type calcium current- CaV1.2),\n"
            "* * iks (IKs current - KCNQ1 + MinK),\n"
            "* * ik1 (IK1 current - KCNN4 a.k.a. KCa3.1),\n"
            "* * ito ([fast] Ito current - Kv4.3 + KChIP2.2).\n"
            "*\n"
            "For each channel you specify drug affinity\n"
            "*   EITHER with IC50 values (in uM), for example for 'herg':\n"
            "* --ic50-herg     hERG   IC50    (optional - defaults to \"no affect\")\n"
            "*   OR with pIC50 values (in log M):\n"
            "* --pic50-herg    hERG   pIC50    (optional - defaults to \"no affect\")\n"
            "*     (you can use a mixture of these for different channels if you wish, \n"
            "*     e.g. --ic50-herg 16600 --pic50-na 5.3 )\n"
            "*   AND specify Hill coefficients (dimensionless):\n"
            "* --hill-herg     hERG   Hill    (optional - defaults to \"1.0\")\n"
            "*\n"
            "* Specifying test concentrations:\n"
            "* --plasma-concs  A list of plasma concentrations at which to test (uM)\n"
            "* OR alternatively:\n"
            "* --plasma-conc-high  Highest plasma concentration to test (uM)\n"
            "* --plasma-conc-low   Lowest  plasma concentration to test (uM) \n"
            "*                                  (optional - defaults to 0)\n"
            "*\n"
            "* both ways of specifying test concentrations have the following optional arguments\n"
            "* --plasma-conc-count  Number of intermediate plasma concentrations to test \n"
            "*                                  (optional - defaults to 0 (for --plasma-concs) or 11 (for --plasma-conc-high)\n"
            "* --plasma-conc-logscale  Use a log spacing for the plasma concentrations \n"
            "*                                  (optional - defaults to false)\n";
    return message;
}

void ApPredictMethods::ReadInIC50AndHill(double& rIc50, double& rHill, const std::string& rChannel)
{
	CommandLineArguments* p_args = CommandLineArguments::Instance();
    if (p_args->OptionExists("--ic50-" + rChannel))
    {
        rIc50 = p_args->GetDoubleCorrespondingToOption("--ic50-" + rChannel);
        if(p_args->OptionExists("--pic50-" + rChannel))
        {
        	EXCEPTION("Duplicate arguments, you cannot specify both IC50 and pIC50 for " << rChannel << " channel.");
        }
    }
    else if(p_args->OptionExists("--pic50-" + rChannel))
    {
        rIc50 = pow(10,-p_args->GetDoubleCorrespondingToOption("--pic50-" + rChannel)+6);
    }

    if (p_args->OptionExists("--hill-" + rChannel))
    {
        rHill = p_args->GetDoubleCorrespondingToOption("--hill-" + rChannel);
    }

    if (!mSuppressOutput) std::cout << "* " << rChannel;
    if (rIc50 > 0)
    {
        if (!mSuppressOutput) std::cout << " IC50 = " << rIc50 << " uM, ";
        if (rHill > 0)
        {
            if (!mSuppressOutput) std::cout << "Hill = " << rHill << "\n";
        }
        else
        {
            if (!mSuppressOutput) std::cout << "Hill = 1.0 (default) \n";
        }
    }
    else
    {
        rHill = -1.0; // Override any values given here to replace with default if there is no channel block.
        if (!mSuppressOutput) std::cout << ": no drug affect\n";
    }
}

void ApPredictMethods::ApplyDrugBlock(boost::shared_ptr<AbstractCvodeCell> pModel,
                                      const std::string& rMetadataName,
                                      const std::string& rShortName,
                                      const double default_conductance,
                                      const double concentration,
                                      const double iC50,
                                      const double hill)
{
    // Here we calculate the proportion of the different channels which are still active
    // (at this concentration of this drug)
    const double conductance_factor = AbstractDataStructure::CalculateConductanceFactor(concentration, iC50,  hill);

    // Some screen output for info.
    if (!mSuppressOutput) std::cout << "g_" << rShortName <<  " factor = " << conductance_factor << "\n";// << std::flush;

    // Check the model has this parameter before we try and set it.
    if (pModel->HasParameter(rMetadataName))
    {
        pModel->SetParameter(rMetadataName, default_conductance*conductance_factor);
    }
    else // We haven't got that conductance parameter, or at least it isn't labelled.
    {
        // If we aren't trying to change it - don't worry, just carry on.
        if (conductance_factor < 1)
        {
            // If the model hasn't got this channel conductance labelled,
            // (but we are trying to change it) throw an error.
            EXCEPTION(pModel->GetSystemName() << " does not have the current \"" << rMetadataName << "\" labelled, but you have requested a block on this channel.");
        }
    }
}

ApPredictMethods::ApPredictMethods()
   : AbstractActionPotentialMethod(),
     mComplete(false)
{
}

void ApPredictMethods::Run()
{
    mProgramName = "Action Potential PreDiCT";
    mOutputFolder = "ApPredict_output/";

    SetupModel setup(this->mHertz); // This class will get model definition from command line.
    mpModel = setup.GetModel();

    CommonRunMethod();
}

void ApPredictMethods::CommonRunMethod()
{
    // Arguments that take default values
    std::vector<double> IC50s;
    std::vector<double> hills;
    std::vector<std::string> metadata_names;
    std::vector<std::string> short_names;

    // Set the initial values of these
    for (unsigned i=0; i<6u; i++)
    {
        IC50s.push_back(-1); // These are our code for 'unset'.
        hills.push_back(-1); // These are our code for 'unset'.
    }

    metadata_names.push_back("membrane_fast_sodium_current_conductance");
    metadata_names.push_back("membrane_L_type_calcium_current_conductance");
    metadata_names.push_back("membrane_rapid_delayed_rectifier_potassium_current_conductance");
    metadata_names.push_back("membrane_slow_delayed_rectifier_potassium_current_conductance");
    metadata_names.push_back("membrane_inward_rectifier_potassium_current_conductance");
    metadata_names.push_back("membrane_fast_transient_outward_current_conductance");

    short_names.push_back("na");
    short_names.push_back("cal");
    short_names.push_back("herg");
    short_names.push_back("iks");
    short_names.push_back("ik1");
    short_names.push_back("ito");

    // Safety Checks
    assert(metadata_names.size()==short_names.size());
    assert(metadata_names.size()==IC50s.size());
    assert(metadata_names.size()==hills.size());

    // Dose calculator asks for some arguments to do with plasma concentrations.
    DoseCalculator dose_calculator;
    mConcs = dose_calculator.GetConcentrations();

    // We check the desired parameters are present in the model, warn if not.
    // This method also changes some metadata names if the model has variants that will do, but aren't ideal
    // and warns if it does this.
    ParameterWrapper(mpModel, metadata_names);

    // Use a helper method to read in IC50 from either --ic50 or --pic50 arguments.
    // Note this is now in micro Molar (1x10^-6 Molar) as per most Pharma use.
    for (unsigned i=0; i<metadata_names.size(); i++)
    {
        ReadInIC50AndHill(IC50s[i], hills[i], short_names[i]);
    }

    if (!mSuppressOutput)
    {
        std::cout << "* max free plasma concentration = " << mConcs.back()<< " uM\n"
                     "* min free plasma concentration = " << mConcs[0] << " uM\n"
                     "* number of plasma concentrations = " << mConcs.size() << "\n";// << std::flush;
    }

    // This just stops the cell automatically resetting the CVODE solver before each call to Solve().
    mpModel->SetAutoReset(false);

    // The following names are fixed and correspond to metadata tags.
    // We record the default parameter values that the model uses.
    // All the drug block models should include these parameter labels
    // But you only get a warning if not, so check the warnings...
    std::vector<double> default_conductances;
    for (unsigned i=0; i<metadata_names.size(); i++)
    {
        double default_value = 1.0;
        if (mpModel->HasParameter(metadata_names[i]))
        {
            default_value = mpModel->GetParameter(metadata_names[i]);
        }
        default_conductances.push_back(default_value);
    }

    boost::shared_ptr<const AbstractOdeSystemInformation> p_ode_info = mpModel->GetSystemInformation();
    std::string model_name = mpModel->GetSystemName();
    boost::static_pointer_cast<RegularStimulus>(mpModel->GetStimulusFunction())->SetStartTime(5.0);

    // Set up foldernames for each model and protocol set.
    std::string foldername = mOutputFolder;

    // Make and clean the above directories.
    mpFileHandler.reset(new OutputFileHandler(foldername));

    // Print out a progress file for monitoring purposes.
    ProgressReporter progress_reporter(foldername, 0.0, (double)(mConcs.size()));
    progress_reporter.PrintInitialising();

    // Open files and write headers
    out_stream steady_voltage_results_file_html = mpFileHandler->OpenOutputFile("voltage_results.html");

    out_stream steady_voltage_results_file = mpFileHandler->OpenOutputFile("voltage_results.dat");
    *steady_voltage_results_file << "Concentration(uM)\tUpstrokeVelocity(mV/ms)\tPeakVm(mV)\tAPD50(ms)\tAPD90(ms)\tdelta_APD90(%)\n";
    *steady_voltage_results_file_html << "<html>\n<head><title>" << mProgramName << " results</title></head>\n";
    *steady_voltage_results_file_html << "<STYLE TYPE=\"text/css\">\n<!--\nTD{font-size: 12px;}\n--->\n</STYLE>\n";
    *steady_voltage_results_file_html << "<body>\n";
    *steady_voltage_results_file_html << "<table width=\"60%\" style=\"background-color:white\" border=\"1\" cellpadding=\"2\" cellspacing=\"0\">\n";
    *steady_voltage_results_file_html << "<tr><td>Concentration (uM)</td><td>Upstroke Velocity (mV/ms)</td><td>Peak membrane voltage (mV)</td><td>APD50 (ms)</td><td>APD90 (ms)</td><td>Change in APD90 (%)</td></tr>\n"; // Header line

    /**
     * START LOOP OVER EACH CONCENTRATION TO TEST WITH
     */
    double control_apd90 = 0;
    for (unsigned conc_index=0; conc_index<mConcs.size(); conc_index++)
    {
        progress_reporter.Update((double)(conc_index));
        std::cout << "Drug Conc = " << mConcs[conc_index] << " uM\n" ;//<< std::flush;

        // Apply drug block on each channel
        for (unsigned i = 0 ; i<metadata_names.size(); i++)
        {
            ApplyDrugBlock(mpModel, metadata_names[i], short_names[i],
                    default_conductances[i], mConcs[conc_index], IC50s[i], hills[i]);
        }

        double apd90, apd50, upstroke, peak;
        OdeSolution solution = SteadyStatePacingExperiment(mpModel, apd90, apd50, upstroke, peak, 0.1 /*ms printing timestep*/, mConcs[conc_index]);

        if (!DidErrorOccur())
        {
            // Record the control APD90 if this concentration is zero.
            if (fabs(mConcs[conc_index])<1e-12)
            {
                control_apd90 = apd90;
            }
            double delta_apd90 = 100*(apd90 - control_apd90)/control_apd90;

            if (!mSuppressOutput) std::cout << mHertz << "Hz Upstroke velocity = " << upstroke << ", Peak mV = " << peak << ", APD50 = " << apd50 << ", APD90 = " << apd90 << ", percent change APD90 = " << delta_apd90 << "\n";// << std::flush;
            *steady_voltage_results_file_html << "<tr><td>"<< mConcs[conc_index] << "</td><td>" << upstroke << "</td><td>" << peak << "</td><td>" << apd50 << "</td><td>" << apd90 << "</td><td>" << delta_apd90 << "</td></tr>\n";
            *steady_voltage_results_file << mConcs[conc_index] << "\t" << upstroke << "\t" << peak << "\t" << apd50  << "\t" << apd90 << "\t" << delta_apd90 << "\n";
        }
        else
        {
            std::string error_code = GetErrorMessage();
            if (!mSuppressOutput) std::cout << mHertz << "Hz Upstroke velocity = " << error_code << ", Peak mV = " << error_code << ", APD50 = " << error_code << ", APD90 = " << error_code << ", percent change APD90 = " << error_code << "\n";// << std::flush;
            *steady_voltage_results_file_html << "<tr><td>"<< mConcs[conc_index] << "</td><td>" << error_code << "</td><td>" << error_code << "</td><td>" << error_code << "</td><td>" << error_code << "</td><td>" << error_code << "</td></tr>\n";
            *steady_voltage_results_file << mConcs[conc_index] << "\t" << error_code << "\t" << error_code << "\t" << error_code  << "\t" << error_code << "\t" << error_code << "\n";
        }

        // Store some things as member variables for returning later (mainly for testing)
        mAPD90s.push_back(apd90); // This is used by TorsadePredict

        // Create unique filename and write the voltage trace to file...
        std::stringstream filename;
        filename << "conc_" << mConcs[conc_index] << "_voltage_trace.dat";
        boost::shared_ptr<RegularStimulus> p_default_stimulus =
                             boost::static_pointer_cast<RegularStimulus>(mpModel->GetStimulusFunction());
        double s1_period = p_default_stimulus->GetPeriod();
        double s_start = p_default_stimulus->GetStartTime();
        std::vector<double> voltages = solution.GetVariableAtIndex(mpModel->GetSystemInformation()->GetStateVariableIndex("membrane_voltage")); // Voltage should always be zero
        ActionPotentialDownsampler(foldername, filename.str(), solution.rGetTimes(), voltages, s1_period, s_start);
    }// Conc

    // Tidy up
    progress_reporter.PrintFinalising();
    *steady_voltage_results_file_html << "</table>\n</body>\n</html>\n";
    steady_voltage_results_file_html->close();
    steady_voltage_results_file->close();

    mComplete = true;
}

void ApPredictMethods::WriteMessageToFile(const std::string& rMessage)
{
    AbstractActionPotentialMethod::WriteMessageToFile(rMessage);
    assert(mpFileHandler);
    // work out what the absolute path of the message file will be
    const std::string message_file_name = mpFileHandler->GetOutputDirectoryFullPath() + "messages.txt";

    // Check to see if this is the first log message
    FileFinder file_finder(message_file_name, RelativeTo::Absolute);
    bool first_message = !file_finder.Exists();

    // Open the messages file, in appending mode.
    out_stream messages_file = mpFileHandler->OpenOutputFile("messages.txt", std::ios::out | std::ios::app);
    if (first_message)
    {
        *messages_file << "Action potential prediction simulation recorded the following notes:\n";
    }
    *messages_file << " * " << rMessage << "\n";
    messages_file->close();
}

std::vector<double> ApPredictMethods::GetConcentrations(void)
{
    if (!mComplete)
    {
        EXCEPTION("Simulation has not been run - check arguments.");
    }
    return mConcs;
}

std::vector<double> ApPredictMethods::GetApd90s(void)
{
    if (!mComplete)
    {
        EXCEPTION("Simulation has not been run - check arguments.");
    }
    return mAPD90s;
}

void ApPredictMethods::ParameterWrapper(boost::shared_ptr<AbstractCvodeCell> pModel, std::vector<std::string>& rMetadataNames)
{
    for (unsigned i=0; i<rMetadataNames.size(); i++)
    {
        if (!pModel->HasParameter(rMetadataNames[i]))
        {
            WARNING(pModel->GetSystemName() << " does not have '" << rMetadataNames[i] << "' labelled, please tag it in the CellML file if it is present.");

            // Not all of the models have a distinct fast I_to component.
            // In this case we look for the complete I_to current instead.
            if (rMetadataNames[i]=="membrane_fast_transient_outward_current_conductance" &&
                pModel->HasParameter("membrane_transient_outward_current_conductance")     )
            {
                WARNING(pModel->GetSystemName() << " does not have 'membrane_fast_transient_outward_current_conductance' labelled, using combined Ito (fast and slow) instead...");
                rMetadataNames[i] = "membrane_transient_outward_current_conductance";
            }
        }
    }
}



