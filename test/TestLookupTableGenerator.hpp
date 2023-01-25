/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTLOOKUPTABLEGENERATOR_HPP_
#define TESTLOOKUPTABLEGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "AbstractUntemplatedLookupTableGenerator.hpp"
#include "LookupTableGenerator.hpp"
#include "NumericFileComparison.hpp"
#include "SetupModel.hpp"
#include "SingleActionPotentialPrediction.hpp"

/**
 * Here we want to generate lookup tables for a given % block of
 * IKr  (hERG)
 * INa  (NaV1.5)
 * ICaL (CaV1.2)
 * IKs  (KCNQ1/minK)
 * Ito  (Kv4.3/KChIP2.2)
 */
class TestLookupTableGenerator : public CxxTest::TestSuite
{
private:
    AbstractUntemplatedLookupTableGenerator *mpGenerator;

public:
    void TestSingleRun()
    {
        SetupModel setup(1.0, 2u); // Ten tusscher '06 at 1 Hz
        boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

        // At this point we could introduce a scaling of the parameters,
        // or change the initial conditions if we have an estimate of the
        // steady state. This is what the lookup table maker will do.

        // Get the whole run packaged into a couple of simple calls
        SingleActionPotentialPrediction ap_prediction(p_model);
        ap_prediction.RunSteadyPacingExperiment();

        // Check it worked
        TS_ASSERT(!ap_prediction.DidErrorOccur());

        // Check some of the results
        TS_ASSERT_DELTA(ap_prediction.GetApd90(), 301.4616, 1e-3);
        TS_ASSERT_DELTA(ap_prediction.GetApd50(), 273.2500, 1e-2);
        TS_ASSERT_DELTA(ap_prediction.GetPeakVoltage(), 37.3489, 1e-2);
        TS_ASSERT_DELTA(ap_prediction.GetUpstrokeVelocity(), 307.6467, 2e-1); // Upstroke sensitive to different versions of CVODE
    }

    void TestLookupTableMaker1d()
    {
        /*
         * For this first test create a 1D hERG block APD90 lookup table.
         */
        unsigned model_index = 2u; // Ten Tusscher 2006 epi

        std::string file_name = "1d_test";
        OutputFileHandler handler("TestLookupTables"); // Wipe the folder for a fresh test each time.

        LookupTableGenerator<1> generator(model_index, file_name, "TestLookupTables");

        TS_ASSERT_THROWS_THIS(generator.SetParameterToScale("sausages", 0.0, 1.0),
                              "tentusscher_model_2006_epi does not have 'sausages' labelled, please tag it in the CellML file.");

        generator.SetParameterToScale("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
        generator.AddQuantityOfInterest(Apd90, 0.1 /*ms*/); // QoI and tolerance

        generator.SetMaxNumEvaluations(5u);
        generator.GenerateLookupTable();

        std::vector<c_vector<double, 1u>> parameter_values = generator.GetParameterPoints();
        std::vector<std::vector<double>> quantities_of_interest = generator.GetFunctionValues();

        TS_ASSERT_EQUALS(parameter_values.size(), 5u);
        TS_ASSERT_EQUALS(quantities_of_interest.size(), 5u);

        for (unsigned i = 0; i < parameter_values.size(); i++)
        {
            std::cout << parameter_values[i][0] << "\t" << quantities_of_interest[i][0] << "\n";
        }

        // Run the generator again
        generator.SetMaxNumEvaluations(10u);
        generator.GenerateLookupTable();

        // Check its new answers
        parameter_values = generator.GetParameterPoints();
        quantities_of_interest = generator.GetFunctionValues();

        TS_ASSERT_EQUALS(parameter_values.size(), 10u);
        TS_ASSERT_EQUALS(quantities_of_interest.size(), 10u);
    }

    void TestLookupTableMaker5d()
    {
        unsigned model_index = 2u; // Ten tusscher '06 (table generated for 1 Hz at present)

        std::string file_name = "5d_test";
        mpGenerator = new LookupTableGenerator<5>(model_index, file_name, "TestLookupTables");

        mpGenerator->SetParameterToScale("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
        mpGenerator->SetParameterToScale("membrane_L_type_calcium_current_conductance", 0.0, 1.0);
        mpGenerator->SetParameterToScale("membrane_fast_sodium_current_conductance", 0.0, 1.0);
        mpGenerator->SetParameterToScale("membrane_slow_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
        mpGenerator->SetParameterToScale("membrane_fast_transient_outward_current_conductance", 0.0, 1.0);

        mpGenerator->AddQuantityOfInterest(Apd90, 0.5 /*ms*/);
        mpGenerator->AddQuantityOfInterest(Apd50, 0.5 /*ms*/);
        mpGenerator->AddQuantityOfInterest(UpstrokeVelocity, 10.0 /* mV/ms */);
        mpGenerator->AddQuantityOfInterest(PeakVoltage, 5 /* mV */);

        mpGenerator->SetMaxNumEvaluations(1u); // This will still do a load of things, for first run, but won't do any refinement.
        mpGenerator->GenerateLookupTable();

        LookupTableGenerator<5> *p_temp = dynamic_cast<LookupTableGenerator<5> *>(mpGenerator);
        std::vector<c_vector<double, 5u>> parameter_values = p_temp->GetParameterPoints();
        std::vector<std::vector<double>> quantities_of_interest = mpGenerator->GetFunctionValues();

        TS_ASSERT_EQUALS(parameter_values.size(), 32u);
        TS_ASSERT_EQUALS(quantities_of_interest.size(), 32u);

        FileFinder human_readable_output("TestLookupTables/" + file_name + ".dat",
                                         RelativeTo::ChasteTestOutput);
        FileFinder human_readable_reference("projects/ApPredict/test/data/" + file_name + ".dat",
                                            RelativeTo::ChasteSourceRoot);
        NumericFileComparison comparer(human_readable_output, human_readable_reference);
        comparer.CompareFiles(1e-4);
    }

    void TestLookupTablesArchiver1d()
    {
        OutputFileHandler handler("TestLookupTableArchiving");
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "Generator1d.arch";

        // Create data structures to store variables to test for equality here
        const unsigned num_evals_before_save = 3u;
        {
            /*
             * For this first test create a 1D hERG block APD90 lookup table.
             */
            unsigned model_index = 2u; // Ten tusscher '06 at 1 Hz
            std::string file_name = "1d_test";

            AbstractUntemplatedLookupTableGenerator *const p_generator = new LookupTableGenerator<1>(model_index, file_name, "TestLookupTableArchiving");

            p_generator->SetParameterToScale("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
            p_generator->AddQuantityOfInterest(Apd90, 0.5 /*ms*/); // QoI and tolerance
            p_generator->SetMaxNumEvaluations(num_evals_before_save);
            p_generator->GenerateLookupTable();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_generator;
            delete p_generator;
        }

        {
            AbstractUntemplatedLookupTableGenerator *p_abstract_generator;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> p_abstract_generator;

            TS_ASSERT_EQUALS(p_abstract_generator->GetDimension(), 1u);

            LookupTableGenerator<1u> *p_generator = dynamic_cast<LookupTableGenerator<1u> *>(p_abstract_generator);

            std::vector<c_vector<double, 1u>> points = p_generator->GetParameterPoints();
            std::vector<std::vector<double>> values = p_generator->GetFunctionValues();

            TS_ASSERT_EQUALS(points.size(), num_evals_before_save);
            TS_ASSERT_EQUALS(values.size(), num_evals_before_save);

            p_generator->SetMaxNumEvaluations(2 * num_evals_before_save);
            p_generator->GenerateLookupTable();

            points = p_generator->GetParameterPoints();
            values = p_generator->GetFunctionValues();

            TS_ASSERT_EQUALS(points.size(), 2 * num_evals_before_save);
            TS_ASSERT_EQUALS(values.size(), 2 * num_evals_before_save);

            delete p_generator;
        }
    }

    void TestLookupTablesArchiver5d()
    {
        OutputFileHandler handler("TestLookupTableArchiving", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "Generator5d.arch";

        // Create data structures to store variables to test for equality here
        const unsigned num_evals_before_save = 32u;
        {
            // Save this generator we have sneakily kept a pointer to from a previous test.
            AbstractUntemplatedLookupTableGenerator *const p_generator = mpGenerator;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_generator;
            delete p_generator; // also deletes mpGenerator...
        }

        {
            AbstractUntemplatedLookupTableGenerator *p_abstract_generator;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> p_abstract_generator;

            TS_ASSERT_EQUALS(p_abstract_generator->GetDimension(), 5u);

            LookupTableGenerator<5u> *p_generator = dynamic_cast<LookupTableGenerator<5u> *>(p_abstract_generator);

            std::vector<c_vector<double, 5u>> points = p_generator->GetParameterPoints();
            std::vector<std::vector<double>> values = p_generator->GetFunctionValues();

            TS_ASSERT_EQUALS(points.size(), num_evals_before_save);
            TS_ASSERT_EQUALS(values.size(), num_evals_before_save);

            delete p_generator;
        }
    }
};

#endif // TESTLOOKUPTABLEGENERATOR_HPP_
