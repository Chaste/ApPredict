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

#ifndef TESTDAVIES2012PAPER_HPP_
#define TESTDAVIES2012PAPER_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCvodeCell.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "CellProperties.hpp"
#include "OutputFileHandler.hpp"
#include "ProgressReporter.hpp"
#include "RegularStimulus.hpp"
#include "SingleActionPotentialPrediction.hpp"
#include "SteadyStateRunner.hpp"

#include "DaviesDogDataStructure.hpp"
#include "davies_isap_2012Cvode.hpp"

class TestDavies2012Paper : public CxxTest::TestSuite
{
public:
    void TestDogDataStructure()
    {
        FileFinder dog_parameter_file("projects/ApPredict/test/data/davies_dog_parameters.txt",
                                      RelativeTo::ChasteSourceRoot);

        DaviesDogDataStructure data(dog_parameter_file);

        TS_ASSERT_EQUALS(data.GetNumDogs(), 20u);

        TS_ASSERT_DELTA(data.GetIKrParameter(0u), 1.00, 1e-12);
        TS_ASSERT_DELTA(data.GetIKrParameter(14u), 3.4373, 1e-4);

        TS_ASSERT_DELTA(data.GetItoGateParameter(19u), 496.1084, 1e-4);
    }

    void TestReproducePaperResults()
    {
        // Set up an isAP model with a CVODE solver.
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
        boost::shared_ptr<AbstractCvodeCell> p_model(new Celldavies_isap_2012FromCellMLCvode(p_solver, p_stimulus));
        boost::shared_ptr<RegularStimulus> p_default_stimulus = p_model->UseCellMLDefaultStimulus();
        p_model->SetTolerances(1e-6, 1e-8); // To be consistent with AZ matlab.

        // Get the different dog parameters out of the data file.
        FileFinder dog_parameter_file("projects/ApPredict/test/data/davies_dog_parameters.txt", RelativeTo::ChasteSourceRoot);
        DaviesDogDataStructure data(dog_parameter_file);

        // The following names are fixed and correspond to metadata tags.
        // We record the default parameter values that the model uses and const them to avoid problems!
        const double default_g_na = p_model->GetParameter("membrane_fast_sodium_current_conductance");
        const double default_g_cal = p_model->GetParameter("membrane_L_type_calcium_current_conductance");
        const double default_g_kr = p_model->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance");
        //const double default_g_ks = p_model->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");
        const double default_g_to = p_model->GetParameter("membrane_transient_outward_current_conductance");
        const double default_g_k1 = p_model->GetParameter("membrane_inward_rectifier_potassium_current_conductance");
        const double default_g_cab = p_model->GetParameter("membrane_background_calcium_current_conductance");
        const double default_g_pca = p_model->GetParameter("membrane_calcium_pump_current_conductance");
        const double default_g_ncx = p_model->GetParameter("membrane_sodium_calcium_exchanger_current_conductance");
        const double default_g_nak = p_model->GetParameter("membrane_sodium_potassium_pump_current_permeability");
        const double default_g_nal = p_model->GetParameter("membrane_persistent_sodium_current_conductance");

        p_default_stimulus->SetStartTime(1.0); // Make the graphs look prettier with this.
        std::string model_name = p_model->GetSystemName();

        unsigned max_num_paces = 4u;
        if (CommandLineArguments::Instance()->OptionExists("--max-num-paces"))
        {
            max_num_paces = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--max-num-paces") - 2u;
        }

        // Set up foldername
        std::string foldername = "TestAstraZenecaModel";

        // Make and clean the above directories.
        OutputFileHandler handler(foldername, false);

        ProgressReporter progress_reporter(foldername, 0.0, (double)(data.GetNumDogs()));
        progress_reporter.PrintInitialising();

        // Open files and write headers
        std::stringstream summary_file_name;
        summary_file_name << "voltage_results_" << max_num_paces + 2 << "_paces.dat";
        out_stream steady_voltage_results_file = handler.OpenOutputFile(summary_file_name.str());
        *steady_voltage_results_file << "Dog\tUpstrokeVelocity(mV/ms)\tPeakVm(mV)\tAPD50(ms)\tAPD90(ms)\tdelta_APD90(%)\n";

        N_Vector solution_at_control = p_model->GetStateVariables(); // Model is at a steady-state for 'base' 1Hz.

        /**
         * START LOOP OVER EACH DOG
         */
        double control_apd90 = 0;
        for (unsigned dog_index = 0; dog_index < data.GetNumDogs(); dog_index++)
        {
            p_model->SetStateVariables(solution_at_control);
            progress_reporter.Update((double)(dog_index));
            std::cout << "DOG = " << dog_index << "\n"
                      << std::flush;

            // Here we calculate the proportion of the different channels which are still active
            // (at this concentration of this drug)
            const double gNa_factor = data.GetINaParameter(dog_index);
            const double gCaL_factor = data.GetICaLParameter(dog_index);
            const double gKr_factor = data.GetIKrParameter(dog_index);
            //const double gKs_factor = data.GetIKsParameter(dog_index);
            const double gto_factor = data.GetItoParameter(dog_index);
            const double gK1_factor = data.GetIK1Parameter(dog_index);
            const double gcab_factor = data.GetICabParameter(dog_index);
            const double pca_factor = data.GetIpCaParameter(dog_index);
            const double ncx_factor = data.GetINcxParameter(dog_index);
            const double nak_factor = data.GetINaKParameter(dog_index);
            const double gnal_factor = data.GetINaLParameter(dog_index);

            std::cout << "gKr factor = " << gKr_factor << "\n";
            std::cout << "gto factor = " << gto_factor << "\n";
            std::cout << "gNa factor = " << gNa_factor << "\n";
            std::cout << "gCaL factor = " << gCaL_factor << "\n";
            //std::cout << "gKs factor = " << gKs_factor << "\n";
            std::cout << "gK1 factor = " << gK1_factor << "\n";
            std::cout << "gcab_factor = " << gcab_factor << "\n";
            std::cout << "pca_factor = " << pca_factor << "\n";
            std::cout << "ncx_factor = " << ncx_factor << "\n";
            std::cout << "nak_factor = " << nak_factor << "\n";
            std::cout << "gnal_factor = " << gnal_factor << "\n";
            std::cout << "d_gate_power_tau = " << data.GetICaLTauPowerParameter(dog_index) << "\n";
            std::cout << "ito_gate_parameter = " << data.GetItoGateParameter(dog_index) << "\n"
                      << std::flush;

            // The following names are fixed and correspond to metadata tags.
            // These ten are all multipliers
            p_model->SetParameter("membrane_fast_sodium_current_conductance", default_g_na * gNa_factor);
            p_model->SetParameter("membrane_L_type_calcium_current_conductance", default_g_cal * gCaL_factor);
            p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", default_g_kr * gKr_factor);
            //p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance",default_g_ks*gKs_factor);
            p_model->SetParameter("membrane_transient_outward_current_conductance", default_g_to * gto_factor);
            p_model->SetParameter("membrane_inward_rectifier_potassium_current_conductance", default_g_k1 * gK1_factor);
            p_model->SetParameter("membrane_background_calcium_current_conductance", default_g_cab * gcab_factor);
            p_model->SetParameter("membrane_calcium_pump_current_conductance", default_g_pca * pca_factor);
            p_model->SetParameter("membrane_sodium_calcium_exchanger_current_conductance", default_g_ncx * ncx_factor);
            p_model->SetParameter("membrane_sodium_potassium_pump_current_permeability", default_g_nak * nak_factor);
            p_model->SetParameter("membrane_persistent_sodium_current_conductance", default_g_nal * gnal_factor);

            // These two are straight parameter swaps (not scalings)
            p_model->SetParameter("membrane_L_type_calcium_current_d_gate_power_tau", data.GetICaLTauPowerParameter(dog_index));
            p_model->SetParameter("membrane_transient_outward_current_time_independent_rectification_gate_constant", data.GetItoGateParameter(dog_index));

            // Use the new class for running to steady state and getting the action potential markers evaluated.
            SingleActionPotentialPrediction ap_prediction(p_model);
            ap_prediction.SetMaxNumPaces(max_num_paces);
            OdeSolution solution = ap_prediction.RunSteadyPacingExperiment();

            if (!ap_prediction.DidErrorOccur())
            {
                double apd90 = ap_prediction.GetApd90();
                double apd50 = ap_prediction.GetApd50();
                double peak = ap_prediction.GetPeakVoltage();
                double upstroke = ap_prediction.GetUpstrokeVelocity();

                // Record the control APD90 if this concentration is zero.
                if (dog_index == 0u)
                {
                    control_apd90 = apd90;
                }
                double delta_apd90 = 100 * (apd90 - control_apd90) / control_apd90;

                std::cout << "Upstroke velocity = " << upstroke << ", Peak mV = " << peak << ", APD50 = " << apd50 << ", APD90 = " << apd90 << ", percent change APD90 = " << delta_apd90 << "\n"; // << std::flush;
                *steady_voltage_results_file << dog_index << "\t" << upstroke << "\t" << peak << "\t" << apd50 << "\t" << apd90 << "\t" << delta_apd90 << "\n";
            }
            else
            {
                std::string error_code = ap_prediction.GetErrorMessage();
                std::cout << "Upstroke velocity = " << error_code << ", Peak mV = " << error_code << ", APD50 = " << error_code << ", APD90 = " << error_code << ", percent change APD90 = " << error_code << "\n"; // << std::flush;
                *steady_voltage_results_file << dog_index << "\t" << error_code << "\t" << error_code << "\t" << error_code << "\t" << error_code << "\t" << error_code << "\n";
            }

            // Create unique filename and write the voltage trace to file...
            std::stringstream filename;
            filename << "dog_" << dog_index << "_voltage_trace_" << max_num_paces + 2u << "_paces.dat";
            {
                out_stream output_file = handler.OpenOutputFile(filename.str());

                *output_file << "Time(ms)\tMembrane_Voltage(mV)\n";
                double s_start = boost::static_pointer_cast<RegularStimulus>(p_model->GetStimulusFunction())->GetStartTime();
                std::vector<double> voltages = solution.GetVariableAtIndex(p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage"));
                double start_time_for_this_pace = solution.rGetTimes()[0] + s_start;
                for (unsigned i = 0; i < voltages.size(); i++)
                {
                    *output_file << solution.rGetTimes()[i] - start_time_for_this_pace << "\t" << voltages[i] << "\n";
                }
                output_file->close();
            }
        } // Dog

        // Tidy up
        progress_reporter.PrintFinalising();
        steady_voltage_results_file->close();
        DeleteVector(solution_at_control);
    }
};

#endif // TESTDAVIES2012PAPER_HPP_
