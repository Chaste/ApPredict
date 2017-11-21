/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef TESTTROUBLESOMEAPEVALUATIONS_HPP_
#define TESTTROUBLESOMEAPEVALUATIONS_HPP_

#include <cxxtest/TestSuite.h>
#include "SetupModel.hpp"
#include "SingleActionPotentialPrediction.hpp"
#include "ZeroStimulus.hpp"

class TestTroublesomeApEvaluations : public CxxTest::TestSuite
{
private:
    std::string Run(boost::shared_ptr<AbstractCvodeCell> pModel,
                    unsigned maxNumPaces = 1800u,
                    bool printTrace = false,
                    double voltageThreshold = DOUBLE_UNSET,
                    double defaultApd90 = DOUBLE_UNSET,
                    std::string fileSuffix = "",
                    double defaultTimePeakVm = DOUBLE_UNSET)
    {
        SingleActionPotentialPrediction runner(pModel);
        runner.SetMaxNumPaces(maxNumPaces);
        runner.SetAlternansIsError();
        runner.SetLackOfOneToOneCorrespondenceIsError();
        if (voltageThreshold != DOUBLE_UNSET)
        {
            runner.SetVoltageThresholdForRecordingAsActionPotential(voltageThreshold);
        }
        if (defaultApd90 != DOUBLE_UNSET)
        {
            runner.SetControlActionPotentialDuration90(defaultApd90);
        }
        if (defaultTimePeakVm != DOUBLE_UNSET)
        {
            runner.SetControlTimeOfPeakVoltage(defaultTimePeakVm);
        }
        OdeSolution soln = runner.RunSteadyPacingExperiment();

        if (printTrace)
        {
            OutputFileHandler handler("Troublesome_AP_debug", false);
            out_stream p_file = handler.OpenOutputFile("voltage_trace" + fileSuffix + ".txt");
            std::vector<double> voltage = soln.GetAnyVariable("membrane_voltage");
            for (unsigned t = 0; t < soln.rGetTimes().size(); t++)
            {
                *p_file << soln.rGetTimes()[t] << "\t" << voltage[t] << std::endl;
            }
            *p_file << std::flush;
            p_file->close();
        }

        if (runner.DidErrorOccur())
        {
            return runner.GetErrorMessage();
        }
        else
        {
            return "No error";
        }
    }

public:
    void TestTroublesomeTenTusscherCases() throw(Exception)
    {
        SetupModel setup(1.0, 2u); // TT06, 1.0Hz
        boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

        const std::string gkr_name = "membrane_rapid_delayed_rectifier_potassium_current_conductance";
        const std::string gks_name = "membrane_slow_delayed_rectifier_potassium_current_conductance";

        const double gKr_max = p_model->GetParameter(gkr_name);
        const double gKs_max = p_model->GetParameter(gks_name);

        std::cout << "Normal model:\n";
        std::string message = Run(p_model);

        TS_ASSERT_EQUALS(message, "No error");

        // We are now in a steady state.
        N_Vector steady_state = p_model->GetStateVariables(); // Take a copy of the state variables.

        std::cout << "Case 1:\n";
        // Alternans, one AP followed by fail to repolarize.
        p_model->SetParameter(gkr_name, gKr_max * 0.0234375);
        p_model->SetParameter(gks_name, gKs_max * 0.046875);
        message = Run(p_model);

        TS_ASSERT_EQUALS(message, "NoActionPotential_2");

        // Alternans, failure to repolarize followed by one AP.
        std::cout << "Case 2:\n";
        p_model->SetStateVariables(steady_state);
        p_model->SetParameter(gkr_name, gKr_max * 0.0390625);
        p_model->SetParameter(gks_name, gKs_max * 0.046875);
        message = Run(p_model);

        TS_ASSERT_EQUALS(message, "NoActionPotential_2");

        // Alternans, two different failures to repolarize.
        std::cout << "Case 3:\n";
        p_model->SetStateVariables(steady_state);
        p_model->SetParameter(gkr_name, gKr_max * 0.265625);
        p_model->SetParameter(gks_name, gKs_max * 0.0234375);
        message = Run(p_model);

        TS_ASSERT_EQUALS(message, "NoActionPotential_2");

        // Lack of 1:1 stimulus APD.
        std::cout << "Case 4:\n";
        p_model->SetStateVariables(steady_state);
        p_model->SetParameter(gkr_name, gKr_max * 0.0);
        p_model->SetParameter(gks_name, gKs_max * 0.075);
        message = Run(p_model);

        TS_ASSERT_EQUALS(message, "NoActionPotential_3");

        // Failure to depolarize (use an effectively ZeroStimulus,
        // very difficult to get TT06 not to fire by blocking currents when threshold is -50 (default)
        std::cout << "Case 5:\n";
        p_model->SetStimulusFunction(boost::shared_ptr<ZeroStimulus>(new ZeroStimulus));
        p_model->SetStateVariables(steady_state);
        p_model->SetParameter(gkr_name, gKr_max);
        p_model->SetParameter(gks_name, gKs_max);

        TS_ASSERT_THROWS_THIS(Run(p_model),
                              "AbstractActionPotentialMethod only works with cells that have a RegularStimulus set.");

        // Make a really weak stimulus.
        p_model->SetStimulusFunction(boost::shared_ptr<RegularStimulus>(new RegularStimulus(-0.01, 3, 1000, 1)));
        message = Run(p_model);

        // So that cell does not depolarize.
        TS_ASSERT_EQUALS(message, "NoActionPotential_1");

        DeleteVector(steady_state);
    }

    void TestTroublesomeOHaraActionPotentials() throw(Exception)
    {
        SetupModel setup(1.0, 6u); // O'Hara, 1.0Hz
        boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

        const std::string gkr_name = "membrane_rapid_delayed_rectifier_potassium_current_conductance";
        const std::string gna_name = "membrane_fast_sodium_current_conductance";

        const double gKr_max = p_model->GetParameter(gkr_name);
        const double gNa_max = p_model->GetParameter(gna_name);
        const double voltage_threshold = -40.42461370307312; // From 10% over max V with gNa=0.

        std::cout << "Normal model:\n";
        std::string message = Run(p_model);
        TS_ASSERT_EQUALS(message, "No error");
        // We are now in a steady state.
        N_Vector steady_state = p_model->GetStateVariables(); // Take a copy of the state variables.

        p_model->SetParameter(gkr_name, gKr_max * 0.5); // For all of these set gKr = 0.5 and vary gNa.

        {
            std::cout << "\nCase 1a: no depolarisation:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gna_name, gNa_max * 0);
            message = Run(p_model, 100, false, voltage_threshold);
            TS_ASSERT_EQUALS(message, "NoActionPotential_1");
        }

        {
            std::cout << "\nCase 1b: no depolarisation again:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gna_name, gNa_max * 0.05);
            message = Run(p_model, 100, false, voltage_threshold);
            TS_ASSERT_EQUALS(message, "NoActionPotential_1");
        }

        {
            std::cout << "\nCase 2: no depolarisation on second AP:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gna_name, gNa_max * 0.065);
            message = Run(p_model, 100, false, voltage_threshold);
            TS_ASSERT_EQUALS(message, "NoActionPotential_5");
        }

        {
            std::cout << "\nCase 3: alternans:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gna_name, gNa_max * 0.081);
            message = Run(p_model, 100, false, voltage_threshold);
            TS_ASSERT_EQUALS(message, "NoActionPotential_4");
        }

        {
            std::cout << "\nCase 4: alternans again 1:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gna_name, gNa_max * 0.092);
            message = Run(p_model, 100, false, voltage_threshold);
            TS_ASSERT_EQUALS(message, "NoActionPotential_4");
        }

        {
            std::cout << "\nCase 5: alternans again 2:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gna_name, gNa_max * 0.093);
            message = Run(p_model, 100, false, voltage_threshold);
            TS_ASSERT_EQUALS(message, "NoActionPotential_4");
        }

        {
            std::cout << "\nCase 6: alternans again 3:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gna_name, gNa_max * 0.094);
            message = Run(p_model, 100, false, voltage_threshold);
            TS_ASSERT_EQUALS(message, "NoActionPotential_4");
        }

        {
            std::cout << "\nCase 7: behaves:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gna_name, gNa_max * 0.1);
            message = Run(p_model, 100, false, voltage_threshold);
            TS_ASSERT_EQUALS(message, "No error");
        }

        const double default_apd = 500.0;

        // Now set sodium and vary gKr in the test cases:
        p_model->SetParameter(gna_name, gNa_max * 0.5);

        {
            std::cout << "\nCase 8a: NoAP3 (but a bit like 5):\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gkr_name, gKr_max * 0.0407);
            message = Run(p_model, 100, true, voltage_threshold, default_apd, "_gKr_0.0407");
            TS_ASSERT_EQUALS(message, "NoActionPotential_3");
        }

        {
            std::cout << "\nCase 8b: NoAP3 (but a bit like 5):\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gkr_name, gKr_max * 0.0408);
            message = Run(p_model, 100, true, voltage_threshold, default_apd, "_gKr_0.0408");
            TS_ASSERT_EQUALS(message, "NoActionPotential_3");
        }

        {
            std::cout << "\nCase 8c: NoAP6 (but a bit like 4 but longer APDs):\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gkr_name, gKr_max * 0.0409);
            message = Run(p_model, 100, true, voltage_threshold, default_apd, "_gKr_0.0409");
            TS_ASSERT_EQUALS(message, "NoActionPotential_6");
        }

        {
            std::cout << "\nCase 9: Depolarisation failure long APDs and alternans:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gkr_name, gKr_max * 0.697631);
            p_model->SetParameter(gna_name, gNa_max * 0.0992804);
            message = Run(p_model, 100, true, voltage_threshold, default_apd, "_gNa_0.0992804_gKr_0.697631");
            // These are long alternans, but caused by depolarisation failure, so should get error code 4.
            TS_ASSERT_EQUALS(message, "NoActionPotential_4");
        }

        const double default_time_of_peak_Vm = 6.0;

        {
            std::cout << "\nCase 10: Only-just-firing sodium channels:\n"
                      << std::endl;
            p_model->SetStateVariables(steady_state);
            p_model->SetParameter(gkr_name, gKr_max * 0.6266);
            p_model->SetParameter(gna_name, gNa_max * 0.1);
            message = Run(p_model, 100, true, voltage_threshold, default_apd, "_gNa_0.1_gKr_0.6266", default_time_of_peak_Vm);
            // These are long alternans, but caused by depolarisation failure, so should get error code 4.
            TS_ASSERT_EQUALS(message, "No error");
        }

        DeleteVector(steady_state);
    }
};

#endif // TESTTROUBLESOMEAPEVALUATIONS_HPP_
