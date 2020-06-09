/*

Copyright (c) 2005-2020, University of Oxford.
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

#include <boost/scoped_array.hpp> // to avoid variable length arrays.
#include <iomanip>                // for setprecision()
#include <pthread.h>              // for pthread_create, pthread_join, etc
#include <unistd.h>               // Timing delays to make pthreads behave themselves

#include "FileFinder.hpp"
#include "LookupTableGenerator.hpp"
#include "ParameterBox.hpp"
#include "SetupModel.hpp"
#include "SingleActionPotentialPrediction.hpp"

void *ThreadedActionPotential(void *argument); // Forward declaration.

struct ThreadReturnData
{
    bool exceptionOccurred;
    std::string exceptionMessage;
    unsigned errorOccurred;
    std::vector<double> QoIs;
};

struct ThreadInputData
{
    std::vector<double> scalings;
    std::vector<std::string> mParameterNames;
    std::vector<double> mUnscaledParameters;
    std::vector<QuantityOfInterest> mQuantitiesToRecord;
    std::vector<double> mInitialConditions;
    unsigned mMaxNumPaces;
    unsigned mModelIndex;
    double mFrequency;
    double mVoltageThreshold;
};

/* Private constructor - just for archiving */
template <unsigned DIM>
LookupTableGenerator<DIM>::LookupTableGenerator()
    : AbstractUntemplatedLookupTableGenerator(),
      mModelIndex(0u),
      mpParentBox(NULL){};

template <unsigned DIM>
LookupTableGenerator<DIM>::LookupTableGenerator(
    const unsigned &rModelIndex, const std::string &rOutputFileName,
    const std::string &rOutputFolder)
    : AbstractUntemplatedLookupTableGenerator(),
      mModelIndex(rModelIndex),
      mFrequency(1.0),
      mMaxNumEvaluations(UNSIGNED_UNSET),
      mNumEvaluations(0u),
      mOutputFileName(rOutputFileName),
      mOutputFolder(rOutputFolder),
      mGenerationHasBegun(false),
      mMaxRefinementDifference(UNSIGNED_UNSET),
      mpParentBox(new ParameterBox<DIM>(NULL)),
      mMaxNumPaces(UNSIGNED_UNSET),
      mVoltageThreshold(-50.0)
{
    // empty
}

template <unsigned DIM>
LookupTableGenerator<DIM>::~LookupTableGenerator()
{
    delete mpParentBox;
}

template <unsigned DIM>
bool LookupTableGenerator<DIM>::GenerateLookupTable()
{
    if (mParameterNames.size() != DIM)
    {
        EXCEPTION(
            "Please add parameter(s) over which to construct a lookup table.");
    }

    if (mQuantitiesToRecord.size() == 0u)
    {
        EXCEPTION(
            "Please add some quantities of interest to construct a lookup table "
            "for.");
    }

    // Get a pointer to the output file for us to use when writing
    OutputFileHandler handler(mOutputFolder, false);
    FileFinder output_file = handler.FindFile(mOutputFileName + ".dat");
    out_stream p_file;

    // Overwrite any existing output file as we will dump stored results from our
    // archive anyway.
    p_file = handler.OpenOutputFile(mOutputFileName + ".dat");

    *p_file << std::setprecision(8);

    // Write out the header line - no longer auto-read, but easy to read by eye so we keep it.
    *p_file << mParameterNames.size() << "\t" << mQuantitiesToRecord.size();
    for (unsigned i = 0; i < mParameterNames.size(); i++)
    {
        *p_file << "\t" << mParameterNames[i];
    }
    for (unsigned i = 0; i < mQuantitiesToRecord.size(); i++)
    {
        // Write out enum as ints
        *p_file << "\t" << (int)(mQuantitiesToRecord[i]);
    }
    *p_file << std::endl;

    // Do a few special things the first time round.
    if (!mGenerationHasBegun)
    {
        std::cout << "Generating from fresh" << std::endl;
        // Provide an initial guess for steady state ICs.
        SetupModel setup(mFrequency, mModelIndex); // model at desired frequency
        boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

        SteadyStateRunner steady_runner(p_model);
        steady_runner.RunToSteadyState();

        // Record these initial conditions (we'll always start from these).
        mInitialConditions = MakeStdVec(p_model->rGetStateVariables());

        // First thing to do is to record the unscaled parameter values.
        for (unsigned i = 0; i < mParameterNames.size(); i++)
        {
            double default_value;
            if (p_model->HasParameter(mParameterNames[i]))
            {
                default_value = p_model->GetParameter(mParameterNames[i]);
            }
            else if (p_model->HasParameter(mParameterNames[i] + "_scaling_factor"))
            {
                default_value = p_model->GetParameter(mParameterNames[i] + "_scaling_factor");
            }
            mUnscaledParameters.push_back(default_value);
        }

        // We now do a special run of a model with sodium current set to zero, so we can see the effect
        // of simply a stimulus current, and then set the threshold for APs accordingly.
        {
            SingleActionPotentialPrediction ap_runner(p_model);
            ap_runner.SuppressOutput();
            ap_runner.SetMaxNumPaces(100u);
            mVoltageThreshold = ap_runner.DetectVoltageThresholdForActionPotential();
        }
        p_model->SetStateVariables(mInitialConditions); // Put the model back to sensible state

        // Initial scalings
        CornerSet set_of_points = mpParentBox->GetCorners();
        assert(set_of_points.size() == pow(2, DIM));

        // Run these initial evaluations multi-threaded.
        RunEvaluationsForThesePoints(set_of_points, p_file);

        mGenerationHasBegun = true;
    }
    else // If generation has already begun then dump the existing results to file.
    // (we are probably recovering an archive and the pre-existing .dat file may be gone).
    {
        std::cout << "Generation has already begun" << std::endl;
        for (unsigned i = 0; i < mParameterPointData.size(); i++)
        {
            std::stringstream line_of_output;
            line_of_output << std::setprecision(8);
            for (unsigned j = 0; j < DIM; j++)
            {
                line_of_output << mParameterPoints[i][j] << "\t";
            }
            line_of_output << mParameterPointData[i]->GetErrorCode();
            for (unsigned j = 0; j < mParameterPointData[i]->rGetQoIs().size(); j++)
            {
                line_of_output << "\t" << mParameterPointData[i]->rGetQoIs()[j];
            }
            if (mParameterPointData[i]->HasErrorEstimates())
            {
                unsigned num_estimates = mParameterPointData[i]->rGetQoIErrorEstimates().size();
                line_of_output << "\t" << num_estimates;
                for (unsigned j = 0; j < num_estimates; j++)
                {
                    line_of_output << "\t"
                                   << mParameterPointData[i]->rGetQoIErrorEstimates()[j];
                }
            }
            *p_file << line_of_output.str() << std::endl;
        }
    }

    bool meets_all_tolerances = false;
    for (unsigned quantitiy_idx = 0u; quantitiy_idx < mQuantitiesToRecord.size();
         quantitiy_idx++)
    {
        // While we are still less than the maximum number of evaluations then
        // refine boxes.
        bool meets_tolerance = false;
        while (mNumEvaluations < mMaxNumEvaluations)
        {
            // Find which parameter box has the largest variation between its corners
            ParameterBox<DIM> *p_box = mpParentBox->FindBoxWithLargestQoIErrorEstimate(
                quantitiy_idx, mQoITolerances[quantitiy_idx],
                mMaxRefinementDifference);

            // If we don't get a box back, then we can quit this while loop,
            // as variation in this QoI is within tols.
            if (!p_box)
            {
                std::cout
                    << "Error estimates are within requested tolerances... finishing\n"
                    << std::flush;
                meets_tolerance = true;
                break;
            }

            // Subdivide this box (NB if we GetCorners() after this,
            // it includes the new points and makes no sense!).
            CornerSet new_parameter_points = p_box->SubDivide();

            // Evaluate at these points.
            RunEvaluationsForThesePoints(new_parameter_points, p_file);
        }

        if (meets_tolerance && quantitiy_idx == 0u)
        {
            meets_all_tolerances = true;
        }

        if (!meets_tolerance)
        {
            meets_all_tolerances = false;
        }
    }

    p_file->close();

    if (meets_all_tolerances)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::RunEvaluationsForThesePoints(
    CornerSet setOfPoints, out_stream &rFile)
{
    // Setup variables to control threading
    unsigned num_threads = setOfPoints.size();
    boost::scoped_array<ThreadInputData> thread_data(new ThreadInputData[num_threads]);
    // struct ThreadInputData thread_data[num_threads];
    boost::scoped_array<void *> answers(new void *[num_threads]);
    int return_code;

    // Generate the threads
    boost::scoped_array<pthread_t> threads(new pthread_t[num_threads]);

    // Create a couple of counters for convenience
    CornerSetIter iter;
    int i;

    /*
     *This loop launches each of the threads.
     */
    for (iter = setOfPoints.begin(), i = 0;
         iter != setOfPoints.end();
         ++iter, ++i)
    {
        std::vector<double> scalings;
        for (unsigned j = 0; j < DIM; j++)
        {
            scalings.push_back((*(*iter))[j]);
        }
        thread_data[i].scalings = scalings;
        thread_data[i].mParameterNames = mParameterNames;
        thread_data[i].mUnscaledParameters = mUnscaledParameters;
        thread_data[i].mQuantitiesToRecord = mQuantitiesToRecord;
        thread_data[i].mInitialConditions = mInitialConditions;
        thread_data[i].mMaxNumPaces = mMaxNumPaces;
        thread_data[i].mModelIndex = mModelIndex;
        thread_data[i].mFrequency = mFrequency;
        thread_data[i].mVoltageThreshold = mVoltageThreshold;

        // std::cout << "Launching pthread[" << i << "]" << std::endl;

        // Launch the ThreadedActionPotential method on this thread
        return_code = pthread_create(&threads[i], NULL, ThreadedActionPotential,
                                     (void *)&thread_data[i]);

        assert(0 == return_code); // Check launch worked OK
        EXCEPT_IF_NOT(0 == return_code);

        // Horrific seg. faults without the below line - bug in p_threads?
        usleep(1e5); // 0.1 second pause to allow thread to launch properly!
    }

    /*
     * This loop gets the answers back from all the threads.
     */
    for (iter = setOfPoints.begin(), i = 0;
         iter != setOfPoints.end();
         ++iter, ++i)
    {
        // Get the answers back
        return_code = pthread_join(threads[i], &answers[i]);
        assert(0 == return_code);
        EXCEPT_IF_NOT(0 == return_code);

        // Translate back from the structs to sensible formats.
        ThreadReturnData *thread_results = (ThreadReturnData *)answers[i];
        if (thread_results->exceptionOccurred)
        {
            EXCEPTION(
                "A thread threw the exception: " << thread_results->exceptionMessage);
        }

        unsigned error_occurred = thread_results->errorOccurred;
        std::vector<double> results = thread_results->QoIs;
        c_vector<double, DIM> *p_scalings = *iter;
        delete thread_results; // Clean up memory

        // Store all the info in the master process and tell boxes about it.
        {
            boost::shared_ptr<ParameterPointData> data = boost::shared_ptr<ParameterPointData>(
                new ParameterPointData(results, error_occurred));

            mParameterPoints.push_back(*p_scalings);
            mParameterPointData.push_back(data);
            mNumEvaluations++;

            // Tell all parameter boxes this information for future refinement.
            mpParentBox->AssignQoIValues(p_scalings, data);
            // This should have updated our error estimates in the ParameterPointData*

            std::stringstream line_of_output;
            line_of_output << std::setprecision(8);
            for (unsigned j = 0; j < DIM; j++)
            {
                line_of_output << (*p_scalings)[j] << "\t";
            }
            line_of_output << error_occurred;
            for (unsigned j = 0; j < results.size(); j++)
            {
                line_of_output << "\t" << results[j];
            }
            if (data->HasErrorEstimates())
            {
                unsigned num_estimates = data->rGetQoIErrorEstimates().size();
                line_of_output << "\t" << num_estimates;
                for (unsigned j = 0; j < num_estimates; j++)
                {
                    line_of_output << "\t" << data->rGetQoIErrorEstimates()[j];
                    if (i == (int)(num_threads - 1u))
                    {
                        // An extra bit of reporting that might be nice can only be called when all boxes have all corner data
                        // i.e. when the last thread has finished, so will appear sporadically in the output!
                        line_of_output << "\t" << mpParentBox->ReportPercentageOfSpaceWhereToleranceIsMetForQoI(mQoITolerances[j], j);
                    }
                }
            }

            *rFile << line_of_output.str() << std::endl;
        }
    }
}

void *ThreadedActionPotential(void *argument)
{
    // bool debugging_on = true;

    struct ThreadInputData *my_data;
    my_data = (struct ThreadInputData *)argument;

    std::vector<double> scalings = my_data->scalings;
    assert(scalings.size() == my_data->mParameterNames.size());

    SetupModel setup(my_data->mFrequency,
                     my_data->mModelIndex); // Ten tusscher '06 at 1 Hz
    boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

    // Do parameter scalings
    for (unsigned i = 0; i < scalings.size(); i++)
    {
        std::string param_name;
        if (p_model->HasParameter(my_data->mParameterNames[i]))
        {
            param_name = my_data->mParameterNames[i];
        }
        else
        {
            param_name = my_data->mParameterNames[i] + "_scaling_factor";
        }
        p_model->SetParameter(param_name,
                              my_data->mUnscaledParameters[i] * (scalings[i]));
    }

    // Reset the state variables to the 'standard' steady state
    N_Vector state_vars = MakeNVector(my_data->mInitialConditions);
    p_model->SetStateVariables(state_vars);

    SingleActionPotentialPrediction ap_runner(p_model);
    ap_runner.SuppressOutput();
    ap_runner.SetMaxNumPaces(my_data->mMaxNumPaces);
    ap_runner.SetLackOfOneToOneCorrespondenceIsError();
    ap_runner.SetVoltageThresholdForRecordingAsActionPotential(
        my_data->mVoltageThreshold);

    // Call the SingleActionPotentialPrediction methods.
    try
    {
        //        if (debugging_on)
        //        {
        //            std::stringstream filename;
        //            for (unsigned i = 0; i < scalings.size(); i++)
        //            {
        //                filename << scalings[i] << "_";
        //            }
        //            OdeSolution solution = ap_runner.RunSteadyPacingExperiment();
        //            solution.WriteToFile("Debugging_Lookup", filename.str(), "ms",
        //            1, false);
        //        }
        //        else
        //        {
        ap_runner.RunSteadyPacingExperiment();
        //        }
    }
    catch (Exception &e)
    {
        ThreadReturnData *return_data = new ThreadReturnData;
        return_data->exceptionOccurred = true;
        return_data->exceptionMessage = e.GetShortMessage();

        DeleteVector(state_vars);
        pthread_exit(return_data);
    }

    unsigned error_occurred = ap_runner.GetErrorCode(); // 0 if there was no error

    // Record the results
    std::vector<double> results;
    for (unsigned i = 0; i < my_data->mQuantitiesToRecord.size(); i++)
    {
        if (ap_runner.DidErrorOccur())
        {
            std::string error_message = ap_runner.GetErrorMessage();
            std::cout << "Lookup table generator reports that " << error_message
                      << "\n"
                      << std::flush;

            // We could use different numerical codes for different errors here if we
            // wanted to, but for QNet all AP errors are just set to -DBL_MAX.
            if (my_data->mQuantitiesToRecord[i] == QNet)
            {
                results.push_back(-DBL_MAX);
                continue;
            }

            // We could use different numerical codes for different errors here if we
            // wanted to.
            if ((error_message == "NoActionPotential_2" || error_message == "NoActionPotential_3") && (my_data->mQuantitiesToRecord[i] == Apd90 || my_data->mQuantitiesToRecord[i] == Apd50))
            {
                // For an APD calculation failure on repolarisation put in the stimulus
                // period.
                double stim_period = boost::static_pointer_cast<RegularStimulus>(
                                         p_model->GetStimulusFunction())
                                         ->GetPeriod();
                results.push_back(stim_period);
            }
            else
            {
                // For everything else (failure to depolarize "NoActionPotential_1")
                // just put in zero for now.
                results.push_back(0.0);
            }
            continue;
        }

        // No error cases
        double temp;
        if (my_data->mQuantitiesToRecord[i] == Apd90)
        {
            temp = ap_runner.GetApd90();
        }
        else if (my_data->mQuantitiesToRecord[i] == Apd50)
        {
            temp = ap_runner.GetApd50();
        }
        else if (my_data->mQuantitiesToRecord[i] == UpstrokeVelocity)
        {
            temp = ap_runner.GetUpstrokeVelocity();
        }
        else if (my_data->mQuantitiesToRecord[i] == PeakVoltage)
        {
            temp = ap_runner.GetPeakVoltage();
        }
        else if (my_data->mQuantitiesToRecord[i] == QNet)
        {
            temp = ap_runner.CalculateQNet();
        }
        results.push_back(temp);
    }

    ThreadReturnData *return_data = new ThreadReturnData;
    return_data->QoIs = results;
    return_data->errorOccurred = error_occurred;
    return_data->exceptionOccurred = false;

    DeleteVector(state_vars);
    pthread_exit(return_data);
}

template <unsigned DIM>
std::vector<c_vector<double, DIM>>
LookupTableGenerator<DIM>::GetParameterPoints()
{
    return mParameterPoints;
}

template <unsigned DIM>
std::vector<std::vector<double>>
LookupTableGenerator<DIM>::GetFunctionValues()
{
    std::vector<std::vector<double>> results;
    for (unsigned i = 0; i < mParameterPointData.size(); i++)
    {
        results.push_back(mParameterPointData[i]->GetQoIs());
    }
    return results;
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::SetParameterToScale(
    const std::string &rMetadataName, const double &rMin, const double &rMax)
{
    if (mParameterNames.size() == DIM)
    {
        EXCEPTION(
            "All parameters have been defined already. You need to expand the "
            "dimension of your Lookup table generator.");
    }

    if (mGenerationHasBegun)
    {
        EXCEPTION(
            "SetParameterToScale cannot be called after GenerateLookupTable.");
    }

    SetupModel setup(1.0, mModelIndex); // model at 1 Hz
    boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

    // The usual case where this is a parameter
    // We'll keep referring to it as a scaling factor, but check before tweaking it whether we need to do this!
    if (p_model->HasParameter(rMetadataName) || p_model->HasParameter(rMetadataName + "_scaling_factor"))
    {
        mParameterNames.push_back(rMetadataName);
    }
    // A special treatment for Ito,fast - we use Ito if it isn't present separately.
    else if (rMetadataName == "membrane_fast_transient_outward_current_conductance" && p_model->HasAnyVariable("membrane_transient_outward_current_conductance"))
    {
        WARNING(p_model->GetSystemName()
                << " does not have "
                   "'membrane_fast_transient_outward_current_conductance' "
                   "labelled, "
                   "using combined Ito (fast and slow) instead...");
        if (p_model->HasParameter("membrane_transient_outward_current_conductance") || p_model->HasParameter("membrane_transient_outward_current_conductance_scaling_factor"))
        {
            mParameterNames.push_back(
                "membrane_transient_outward_current_conductance");
        }
    }
    else // It's neither named nor a scaling factor.
    {
        EXCEPTION(p_model->GetSystemName()
                  << " does not have '" << rMetadataName
                  << "' labelled, please tag it in the CellML file.");
    }

    mMinimums.push_back(rMin);
    mMaximums.push_back(rMax);
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::AddQuantityOfInterest(
    QuantityOfInterest quantity, double tolerance)
{
    if (mGenerationHasBegun)
    {
        EXCEPTION(
            "AddQuantityOfInterest cannot be called after GenerateLookupTable.");
    }
    mQuantitiesToRecord.push_back(quantity);
    mQoITolerances.push_back(tolerance);
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::SetMaxNumEvaluations(
    const unsigned &rMaxNumEvals)
{
    mMaxNumEvaluations = rMaxNumEvals;
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::SetMaxVariationInRefinement(
    const unsigned &rMaxRefinementDifference)
{
    mMaxRefinementDifference = rMaxRefinementDifference;
}

template <unsigned DIM>
std::vector<std::vector<double>> LookupTableGenerator<DIM>::Interpolate(
    const std::vector<c_vector<double, DIM>> &rParameterPoints)
{
    std::vector<std::vector<double>> interpolated_values;

    for (unsigned i = 0; i < rParameterPoints.size(); i++)
    {
        interpolated_values.push_back(
            mpParentBox->InterpolateQoIsAt(rParameterPoints[i]));
    }

    return interpolated_values;
}

template <unsigned DIM>
std::vector<std::vector<double>> LookupTableGenerator<DIM>::Interpolate(
    const std::vector<std::vector<double>> &rParameterPoints)
{
    // Convert std::vector to c_vector...
    std::vector<c_vector<double, DIM>> c_vec_parameter_points;
    for (unsigned i = 0; i < rParameterPoints.size(); i++)
    {
        assert(rParameterPoints[i].size() == DIM);
        c_vector<double, DIM> c_vec_point;
        for (unsigned j = 0; j < DIM; j++)
        {
            c_vec_point[j] = rParameterPoints[i][j];
        }
        c_vec_parameter_points.push_back(c_vec_point);
    }
    // Now just call the method above.
    return Interpolate(c_vec_parameter_points);
}

template <unsigned DIM>
unsigned LookupTableGenerator<DIM>::GetNumEvaluations()
{
    return mNumEvaluations;
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::SetMaxNumPaces(unsigned numPaces)
{
    mMaxNumPaces = numPaces;
}

template <unsigned DIM>
unsigned LookupTableGenerator<DIM>::GetMaxNumPaces()
{
    return mMaxNumPaces;
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::SetPacingFrequency(double frequency)
{
    mFrequency = frequency;
}

template <unsigned DIM>
unsigned LookupTableGenerator<DIM>::GetDimension() const
{
    return DIM;
}

template <unsigned DIM>
std::vector<std::string> LookupTableGenerator<DIM>::GetParameterNames() const
{
    assert(mParameterNames.size() == DIM);
    return mParameterNames;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////
template class LookupTableGenerator<1u>;
template class LookupTableGenerator<2u>;
template class LookupTableGenerator<3u>;
template class LookupTableGenerator<4u>;
template class LookupTableGenerator<5u>;
template class LookupTableGenerator<6u>;
template class LookupTableGenerator<7u>;
// Just up to 7D for now, may need to be bigger eventually.

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS1(LookupTableGenerator, 1u)
EXPORT_TEMPLATE_CLASS1(LookupTableGenerator, 2u)
EXPORT_TEMPLATE_CLASS1(LookupTableGenerator, 3u)
EXPORT_TEMPLATE_CLASS1(LookupTableGenerator, 4u)
EXPORT_TEMPLATE_CLASS1(LookupTableGenerator, 5u)
EXPORT_TEMPLATE_CLASS1(LookupTableGenerator, 6u)
EXPORT_TEMPLATE_CLASS1(LookupTableGenerator, 7u)
