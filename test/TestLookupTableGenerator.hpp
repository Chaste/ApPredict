/*

Copyright (c) 2005-2026, University of Oxford.
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

#include <atomic>
#include <chrono>
#include <cmath>
#include <set>
#include <thread>

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

    void TestLookupTableMaker2dBisection()
    {
        unsigned model_index = 2u; // Ten Tusscher 2006 epi

        std::string file_name = "2d_test";
        LookupTableGenerator<2> generator(model_index, file_name, "TestLookupTables");
        generator.SetParameterToScale("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
        generator.SetParameterToScale("membrane_L_type_calcium_current_conductance", 0.0, 1.0);
        generator.AddQuantityOfInterest(Apd90, 0.5 /*ms*/); // QoI and tolerance

        generator.SetMaxNumEvaluations(1u); // Just does the corners
        generator.GenerateLookupTable();
        TS_ASSERT_EQUALS(generator.GetNumEvaluations(), 4u);

        // Each refinement step bisects one box, adding at most 2^(DIM-1) = 2 new points.
        for (unsigned i = 0; i < 6u; i++)
        {
            unsigned num_evals_before = generator.GetNumEvaluations();
            generator.SetMaxNumEvaluations(num_evals_before + 1u);
            generator.GenerateLookupTable();
            TS_ASSERT_LESS_THAN_EQUALS(num_evals_before + 1u, generator.GetNumEvaluations());
            TS_ASSERT_LESS_THAN_EQUALS(generator.GetNumEvaluations(), num_evals_before + 2u);
        }

        std::vector<c_vector<double, 2u>> parameter_values = generator.GetParameterPoints();
        std::vector<std::vector<double>> quantities_of_interest = generator.GetFunctionValues();
        TS_ASSERT_EQUALS(parameter_values.size(), generator.GetNumEvaluations());
        TS_ASSERT_EQUALS(quantities_of_interest.size(), generator.GetNumEvaluations());

        // Interpolation at the corners of parameter space (the first points evaluated) is exact,
        // and bilinear interpolation shouldn't go outside the range of the data anywhere.
        std::vector<std::vector<double>> interpolated = generator.Interpolate(parameter_values);
        double min_apd = DBL_MAX;
        double max_apd = -DBL_MAX;
        for (unsigned i = 0; i < quantities_of_interest.size(); i++)
        {
            min_apd = std::min(min_apd, quantities_of_interest[i][0]);
            max_apd = std::max(max_apd, quantities_of_interest[i][0]);
            if (i < 4u)
            {
                TS_ASSERT_DELTA(interpolated[i][0], quantities_of_interest[i][0], 1e-12);
            }
        }
        std::vector<std::vector<double>> sample_points;
        for (unsigned i = 0; i <= 10u; i++)
        {
            sample_points.push_back(std::vector<double>{ 0.1 * i, 1.0 - 0.1 * i });
            sample_points.push_back(std::vector<double>{ 0.1 * i, 0.33 });
        }
        interpolated = generator.Interpolate(sample_points);
        for (unsigned i = 0; i < interpolated.size(); i++)
        {
            TS_ASSERT_LESS_THAN_EQUALS(min_apd - 1e-9, interpolated[i][0]);
            TS_ASSERT_LESS_THAN_EQUALS(interpolated[i][0], max_apd + 1e-9);
        }

        // Check archiving and resuming work on a bisected table.
        OutputFileHandler handler("TestLookupTableArchiving", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "Generator2d.arch";
        {
            AbstractUntemplatedLookupTableGenerator* const p_generator = &generator;
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_generator;
        }
        {
            AbstractUntemplatedLookupTableGenerator* p_generator;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_generator;

            TS_ASSERT_EQUALS(p_generator->GetDimension(), 2u);
            TS_ASSERT_EQUALS(p_generator->GetNumEvaluations(), generator.GetNumEvaluations());
            std::vector<std::vector<double>> loaded_interpolated = p_generator->Interpolate(sample_points);
            for (unsigned i = 0; i < interpolated.size(); i++)
            {
                TS_ASSERT_DELTA(loaded_interpolated[i][0], interpolated[i][0], 1e-12);
            }

            unsigned num_evals_before = p_generator->GetNumEvaluations();
            p_generator->SetMaxNumEvaluations(num_evals_before + 1u);
            p_generator->GenerateLookupTable();
            TS_ASSERT_LESS_THAN_EQUALS(num_evals_before + 1u, p_generator->GetNumEvaluations());
            TS_ASSERT_LESS_THAN_EQUALS(p_generator->GetNumEvaluations(), num_evals_before + 2u);
            delete p_generator;
        }
    }

    /**
     * Refine a 2D table of a cheap analytic function, which has a region where 'evaluation' reports
     * an error code, and takes a pseudo-random time to evaluate so that evaluations finish
     * out of order.
     *
     * The function value is continuous across the error region, so boxes on its edge can converge
     * (a jump in value there would be refined down to the minimum box width).
     *
     * @param rGenerator  the generator to set up (parameters must be set already)
     * @param rNumCalls  incremented on every call to the function (from any thread)
     */
    void SetAnalyticEvaluationFunction(LookupTableGenerator<2>& rGenerator, std::atomic<unsigned>& rNumCalls)
    {
        rGenerator.mEvaluationFunctionForTesting = [&rNumCalls](const std::vector<double>& rX, std::vector<double>& rQoIs, unsigned& rErrorCode)
        {
            rNumCalls++;
            // A delay of 0-9ms that varies between points.
            unsigned delay_ms = (unsigned)(std::fabs(std::sin(1000.0 * rX[0] + 37.0 * rX[1])) * 10.0);
            std::this_thread::sleep_for(std::chrono::milliseconds(delay_ms));

            rQoIs.clear();
            rQoIs.push_back(300.0 + 100.0 / (0.2 + rX[0] + 0.5 * rX[1]));
            rErrorCode = (rX[0] + rX[1] < 0.3) ? 2u : 0u;
        };
    }

    void TestParallelRefinementWithAnalyticFunction()
    {
        std::atomic<unsigned> num_calls(0u);
        LookupTableGenerator<2> generator(2u, "2d_analytic_parallel", "TestLookupTableParallel");
        generator.SetParameterToScale("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
        generator.SetParameterToScale("membrane_slow_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
        generator.AddQuantityOfInterest(Apd90, 2.0 /*ms*/);
        generator.SetMaxVariationInRefinement(3u);
        TS_ASSERT_THROWS_THIS(generator.SetNumThreads(0u), "The number of threads must be at least one.");
        generator.SetNumThreads(8u);
        TS_ASSERT_EQUALS(generator.GetNumThreads(), 8u);
        SetAnalyticEvaluationFunction(generator, num_calls);

        // Stop part way through, the cap can only be exceeded by the points on one new plane (1 in 2D).
        generator.SetMaxNumEvaluations(100u);
        TS_ASSERT_EQUALS(generator.GenerateLookupTable(), false);
        TS_ASSERT_LESS_THAN_EQUALS(100u, generator.GetNumEvaluations());
        TS_ASSERT_LESS_THAN_EQUALS(generator.GetNumEvaluations(), 101u);
        CheckGeneratorIsConsistent(generator, num_calls);

        // Archive it, and carry on to convergence with a different number of threads.
        OutputFileHandler handler("TestLookupTableParallel", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "Generator2dParallel.arch";
        {
            AbstractUntemplatedLookupTableGenerator* const p_generator = &generator;
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_generator;
        }
        AbstractUntemplatedLookupTableGenerator* p_abstract_generator;
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_abstract_generator;
        }
        LookupTableGenerator<2>* p_loaded = dynamic_cast<LookupTableGenerator<2>*>(p_abstract_generator);
        TS_ASSERT_EQUALS(p_loaded->GetNumEvaluations(), generator.GetNumEvaluations());
        p_loaded->SetNumThreads(3u);
        SetAnalyticEvaluationFunction(*p_loaded, num_calls);
        p_loaded->SetMaxNumEvaluations(100000u);
        TS_ASSERT_EQUALS(p_loaded->GenerateLookupTable(), true);
        std::cout << "Analytic 2D table converged with " << p_loaded->GetNumEvaluations() << " evaluations.\n";
        TS_ASSERT_LESS_THAN(101u, p_loaded->GetNumEvaluations());
        CheckGeneratorIsConsistent(*p_loaded, num_calls);

        // All the boxes meet the tolerance.
        std::vector<ParameterBox<2>*> leaves;
        p_loaded->CollectLeafBoxes(p_loaded->mpParentBox, leaves);
        for (unsigned i = 0; i < leaves.size(); i++)
        {
            TS_ASSERT_EQUALS(leaves[i]->DoesBoxNeedFurtherRefinement(2.0, 0u), false);
        }

        // And the interpolation is good (looking away from the steepest part of the function).
        for (unsigned i = 0; i <= 20u; i++)
        {
            for (unsigned j = 0; j <= 20u; j++)
            {
                std::vector<double> x{ 0.05 * i, 0.05 * j };
                if (x[0] + x[1] < 0.45)
                {
                    continue;
                }
                double interpolated = p_loaded->Interpolate(std::vector<std::vector<double> >{ x })[0][0];
                TS_ASSERT_DELTA(interpolated, 300.0 + 100.0 / (0.2 + x[0] + 0.5 * x[1]), 4.0);
            }
        }
        delete p_loaded;
    }

    void TestSingleThreadedRefinementIsDeterministic()
    {
        std::vector<std::vector<c_vector<double, 2u> > > runs;
        for (unsigned run = 0; run < 2u; run++)
        {
            std::atomic<unsigned> num_calls(0u);
            LookupTableGenerator<2> generator(2u, "2d_analytic_serial", "TestLookupTableParallel");
            generator.SetParameterToScale("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
            generator.SetParameterToScale("membrane_slow_delayed_rectifier_potassium_current_conductance", 0.0, 1.0);
            generator.AddQuantityOfInterest(Apd90, 2.0 /*ms*/);
            generator.SetNumThreads(1u);
            SetAnalyticEvaluationFunction(generator, num_calls);
            generator.SetMaxNumEvaluations(150u);
            generator.GenerateLookupTable();
            CheckGeneratorIsConsistent(generator, num_calls);
            runs.push_back(generator.GetParameterPoints());
        }
        TS_ASSERT_EQUALS(runs[0].size(), runs[1].size());
        for (unsigned i = 0; i < std::min(runs[0].size(), runs[1].size()); i++)
        {
            TS_ASSERT_DELTA(runs[0][i][0], runs[1][i][0], 1e-12);
            TS_ASSERT_DELTA(runs[0][i][1], runs[1][i][1], 1e-12);
        }
    }

    /**
     * Check that every point was evaluated exactly once, and every box has all its corners evaluated.
     */
    void CheckGeneratorIsConsistent(LookupTableGenerator<2>& rGenerator, std::atomic<unsigned>& rNumCalls)
    {
        std::vector<c_vector<double, 2u> > points = rGenerator.GetParameterPoints();
        TS_ASSERT_EQUALS(points.size(), rGenerator.GetNumEvaluations());
        TS_ASSERT_EQUALS(rNumCalls.load(), rGenerator.GetNumEvaluations());
        std::set<c_vector<double, 2u>*, c_vector_compare<2u> > unique_points;
        for (unsigned i = 0; i < points.size(); i++)
        {
            unique_points.insert(&points[i]);
        }
        TS_ASSERT_EQUALS(unique_points.size(), points.size());

        std::vector<ParameterBox<2>*> leaves;
        rGenerator.CollectLeafBoxes(rGenerator.mpParentBox, leaves);
        for (unsigned i = 0; i < leaves.size(); i++)
        {
            TS_ASSERT(leaves[i]->mAllCornersEvaluated);
            TS_ASSERT_EQUALS(leaves[i]->mParameterPointDataMapPredictions.size(), 0u);
        }
        TS_ASSERT_EQUALS(rGenerator.mpParentBox->GetCorners().size(), rGenerator.GetNumEvaluations());
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
