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
#ifdef CHASTE_CVODE

#ifndef _TESTAPPREDICT_HPP_
#define _TESTAPPREDICT_HPP_

#include <boost/assign/list_of.hpp>
#include <cxxtest/TestSuite.h>

#include <boost/shared_ptr.hpp>
#include "ApPredictMethods.hpp"
#include "CommandLineArgumentsMocker.hpp"
#include "FileFinder.hpp"
#include "NumericFileComparison.hpp"
#include "SetupModel.hpp"

class TestApPredict : public CxxTest::TestSuite
{
public:
    /**
     *
     * This test will wipe $CHASTE_TEST_OUTPUT/ApPredict_output/
     *
     * The first test overwrites CommandLineArguments and checks exceptions are thrown correctly.
     */
    void TestSomeExceptions(void)
    {
        // Check some exceptions are thrown correctly...
        // N.B. the constructor does some of the argument reading, so that needs
        // to be called after the arugment wrapper.
        {
            CommandLineArgumentsMocker wrapper("--plasma-concs 1 10 --pic50-herg 3");

            ApPredictMethods methods;

            TS_ASSERT_THROWS_THIS(methods.Run(),
                                  "Argument \"--model <index or name or file>\" is required (run ApPredict executable with no options for help message).");
        }

        {
            CommandLineArgumentsMocker wrapper("--model 2");

            ApPredictMethods methods;

            TS_ASSERT_THROWS_THIS(methods.Run(),
                                  "Argument \"--plasma-conc-high <concentration in uM>\" or \"--plasma-concs <concentrations in uM>\" is required");
        }

        {
            CommandLineArgumentsMocker wrapper("--model 1 --pacing-freq 0 --pacing-max-time 20 --plasma-concs 1 ");

            TS_ASSERT_THROWS_THIS(ApPredictMethods methods,
                                  "The pacing frequency (0) set by '--pacing-freq' option must be a positive number.");
        }

        {
            CommandLineArgumentsMocker wrapper("--model 1 --cellml 1 --pacing-freq 1 --pacing-max-time 20 --plasma-concs 1 ");

            TS_ASSERT_THROWS_THIS(SetupModel setup(1.0, UNSIGNED_UNSET),
                                  "You can only call ApPredict with the option '--model' OR '--cellml <file>' (not both).");
        }

        {
            CommandLineArgumentsMocker wrapper("--cellml 1 --pacing-freq 1 --pacing-max-time 20 --plasma-concs 1 ");

            TS_ASSERT_THROWS_THIS(SetupModel setup(1.0, UNSIGNED_UNSET),
                                  "Invalid file given with --cellml argument: 1");
        }
        {
            CommandLineArgumentsMocker wrapper("--model bla --pacing-freq 1 --pacing-max-time 20 --plasma-concs 1 ");

            TS_ASSERT_THROWS_THIS(SetupModel setup(1.0, UNSIGNED_UNSET),
                                  "No model matches this index: bla");
        }
        {
            CommandLineArgumentsMocker wrapper("--model 99999");

            TS_ASSERT_THROWS_THIS(SetupModel setup(1.0, UNSIGNED_UNSET),
                                  "No model matches this index: 99999");
        }
        {
            CommandLineArgumentsMocker wrapper("--cellml projects/ApPredict/src/cellml/cellml/ten_tusscher_model_2006_epi.cellml --plasma-concs 1 10 --pic50-herg 4.5 --plasma-conc-logscale false --output-dir ApPredict_output_long");

            ApPredictMethods methods;
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Argument --cellml <file> is deprecated: use --model <file> instead.");
        }
        {
            CommandLineArgumentsMocker wrapper("--cellml bla.cellml");

            TS_ASSERT_THROWS_THIS(SetupModel setup(1.0, UNSIGNED_UNSET),
                                  "Invalid file given with --cellml argument: bla.cellml");
        }
    }

    void TestVoltageThresholdDetectionAlgorithm()
    {
        std::vector<double> thresholds_for_each_model = boost::assign::list_of(-46.7750) /*Shannon etc.*/
            (-23.0772)(-34.6525)(-35.9230)(-28.4091)(-38.4384)(-40.6058);

        for (unsigned model_index = 1; model_index < 8u; model_index++)
        {
            SetupModel setup(1.0, model_index); // models at 1 Hz
            boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

            SingleActionPotentialPrediction ap_runner(p_model);
            ap_runner.SuppressOutput();
            ap_runner.SetMaxNumPaces(100u);
            double threshold_voltage = ap_runner.DetectVoltageThresholdForActionPotential();

            TS_ASSERT_DELTA(threshold_voltage, thresholds_for_each_model[model_index - 1u], 1e-2);
        }
    }

    /**
     * This test should emulate the standalone executable and read your command line arguments.
     */
    void TestDrugAffectByVaryingConductances(void)
    {
        //////////// DEFINE PARAMETERS ///////////////
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments, and
        // char **argv = *(p_args->p_argv); // is a char** of them.
        unsigned num_args = argc - 1;
        std::cout << "# " << num_args << " arguments supplied.\n"
                  << std::flush;

        if (num_args == 0 || CommandLineArguments::Instance()->OptionExists("--help"))
        {
            std::cerr << ApPredictMethods::PrintArguments() << std::flush;
            return;
        }
        ApPredictMethods methods; // No Torsade predictions.
        methods.Run();
    }

    void TestChangingSimulusDuration(void)
    {
        {
            CommandLineArgumentsMocker wrapper("--model 4 --pacing-freq 1 --plasma-concs 0 --pacing-max-time 0.2 --no-downsampling");

            ApPredictMethods methods;

            methods.Run();

            FileFinder generated_file("ApPredict_output/conc_0_voltage_trace.dat", RelativeTo::ChasteTestOutput);
            FileFinder reference_file("projects/ApPredict/test/data/hund_rudy_default_stimulus.dat", RelativeTo::ChasteSourceRoot);
            TS_ASSERT(generated_file.IsFile());
            TS_ASSERT(reference_file.IsFile());

            NumericFileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles(1.5e-2));

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(), 2u);
            TS_ASSERT_DELTA(apd90s[0], 232.777, 2e-2);
            TS_ASSERT_DELTA(apd90s[1], 232.744, 2e-2);
        }
        {
            CommandLineArgumentsMocker wrapper("--model 4 --pacing-freq 1 --plasma-concs 0 --pacing-max-time 0.2 --no-downsampling --pacing-stim-duration 5 --pacing-stim-magnitude -16");

            ApPredictMethods methods;
            methods.SetOutputDirectory("ApPredict_output2/");
            methods.Run();

            FileFinder generated_file("ApPredict_output2/conc_0_voltage_trace.dat", RelativeTo::ChasteTestOutput);
            FileFinder reference_file("projects/ApPredict/test/data/hund_rudy_modified_stimulus.dat", RelativeTo::ChasteSourceRoot);
            TS_ASSERT(generated_file.IsFile());
            TS_ASSERT(reference_file.IsFile());

            NumericFileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles(1.5e-2));

            std::vector<double> apd90s = methods.GetApd90s();
            TS_ASSERT_EQUALS(apd90s.size(), 2u);
            TS_ASSERT_DELTA(apd90s[0], 198.381, 2e-2);
            TS_ASSERT_DELTA(apd90s[1], 198.3093, 2e-2);
        }
    }

    void TestCrash(void)
    {
        CommandLineArgumentsMocker wrapper("--pic50-herg 6 --pic50-spread-herg 0.2 --plasma-concs 10 --credible-intervals --model 8 --pacing-freq 1 --pacing-max-time 5");

        ApPredictMethods methods;

        methods.Run();
    }
};

#endif //_TESTAPPREDICT_HPP_

#endif //_CHASTE_CVODE
