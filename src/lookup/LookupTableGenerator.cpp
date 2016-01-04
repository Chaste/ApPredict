/*

Copyright (c) 2005-2016, University of Oxford.
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

#include <iomanip> // for setprecision()
#include <boost/scoped_array.hpp> // to avoid variable length arrays.

#include "LookupTableGenerator.hpp"
#include "SingleActionPotentialPrediction.hpp"
#include "FileFinder.hpp"
#include "SetupModel.hpp"
#include "ParameterBox.hpp"

void* ThreadedActionPotential(void *argument); // Forward declaration.

struct ThreadReturnData
{
    bool exceptionOccurred;
    std::string exceptionMessage;
    bool errorOccurred;
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
};

template<unsigned DIM>
LookupTableGenerator<DIM>::LookupTableGenerator(const unsigned& rModelIndex,
                                                const std::string& rOutputFileName,
                                                const std::string& rOutputFolder)
    : mModelIndex(rModelIndex),
      mFrequency(1.0),
      mMaxNumEvaluations(UNSIGNED_UNSET),
      mNumEvaluations(0u),
      //mpSingleApRunner(NULL),
      mOutputFileName(rOutputFileName),
      mOutputFolder(rOutputFolder),
      mGenerationHasBegun(false),
      mMaxRefinementDifference(UNSIGNED_UNSET),
      mpParentBox(new ParameterBox<DIM>(NULL)),
      mMaxNumPaces(UNSIGNED_UNSET)
{
    // empty
}

template<unsigned DIM>
LookupTableGenerator<DIM>::~LookupTableGenerator()
{
    delete mpParentBox;
}

template<unsigned DIM>
void LookupTableGenerator<DIM>::GenerateLookupTable()
{
    if (mParameterNames.size() != DIM)
    {
        EXCEPTION("Please add parameter(s) over which to construct a lookup table.");
    }

    if (mQuantitiesToRecord.size() == 0u)
    {
        EXCEPTION("Please add some quantities of interest to construct a lookup table for.");
    }

    // Get a pointer to the output file for us to use when writing
    OutputFileHandler handler(mOutputFolder, false);
    FileFinder output_file = handler.FindFile(mOutputFileName + ".dat");
    out_stream p_file;

    // Overwrite any existing output file as we will dump stored results from our archive anyway.
    p_file = handler.OpenOutputFile(mOutputFileName + ".dat");

    *p_file << std::setprecision(8);

    // Write out the header line - understood by LookupTableReader
    *p_file << mParameterNames.size() << "\t" << mQuantitiesToRecord.size() << "\t";
    for (unsigned i=0; i<mParameterNames.size(); i++)
    {
        *p_file << mParameterNames[i] << "\t";
    }
    for (unsigned i=0; i<mQuantitiesToRecord.size(); i++)
    {
        // Write out enum as ints and then convert back in LookupTableReader.
        *p_file << (int)(mQuantitiesToRecord[i]) << "\t";
    }
    *p_file << std::endl;

    // Do a few special things the first time round.
    if (!mGenerationHasBegun)
    {
        // Provide an initial guess for steady state ICs.
        SetupModel setup(1.0, mModelIndex); // model at 1 Hz
        boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

        SteadyStateRunner steady_runner(p_model);
        steady_runner.RunToSteadyState();

        // Record these initial conditions (we'll always start from these).
        mInitialConditions = MakeStdVec(p_model->rGetStateVariables());

        // First thing to do is to record the unscaled parameter values.
        for (unsigned i=0; i<mParameterNames.size(); i++)
        {
            mUnscaledParameters.push_back(p_model->GetParameter(mParameterNames[i]));
        }

        // Initial scalings
        CornerSet set_of_points = mpParentBox->GetCorners();
        assert(set_of_points.size()==pow(2,DIM));

        // Run these initial evaluations multi-threaded.
        RunEvaluationsForThesePoints(set_of_points, p_file);

        mGenerationHasBegun = true;
    }
    else // If generation has already begun then dump the existing results to file.
         // (we are probably recovering an archive and the pre-existing .dat file may be gone.
    {
        for (unsigned i=0; i<mParameterPointData.size(); i++)
        {
            std::stringstream line_of_output;
            line_of_output << std::setprecision(8);
            for (unsigned j=0; j<DIM; j++)
            {
                line_of_output << mParameterPoints[i][j] << "\t";
            }
            line_of_output << mParameterPointData[i]->DidErrorOccur() << "\t";
            for (unsigned j=0; j< mParameterPointData[i]->rGetQoIs().size(); j++)
            {
                line_of_output << mParameterPointData[i]->rGetQoIs()[j] << "\t";
            }
            if (mParameterPointData[i]->HasErrorEstimates())
            {
                unsigned num_estimates = mParameterPointData[i]->rGetQoIErrorEstimates().size();
                line_of_output << num_estimates << "\t";
                for (unsigned j=0; j< num_estimates; j++)
                {
                    line_of_output << mParameterPointData[i]->rGetQoIErrorEstimates()[j] << "\t";
                }
            }
            *p_file << line_of_output.str() << std::endl;
        }
    }

    for (unsigned quantitiy_idx = 0u; quantitiy_idx<mQuantitiesToRecord.size(); quantitiy_idx++)
    {
        // While we are still less than the maximum number of evaluations then
        // refine boxes.
        while (mNumEvaluations < mMaxNumEvaluations)
        {
            // Find which parameter box has the largest variation between its corners
            ParameterBox<DIM>* p_box =
                    mpParentBox->FindBoxWithLargestQoIErrorEstimate(quantitiy_idx, mQoITolerances[quantitiy_idx], mMaxRefinementDifference);

            // If we don't get a box back, then we can quit this while loop,
            // as variation in this QoI is within tols.
            if (!p_box)
            {
                std::cout << "Error estimates are within requested tolerances... finishing\n" << std::flush;
                break;
            }

            // Subdivide this box (NB if we GetCorners() after this,
            // it includes the new points and makes no sense!).
            CornerSet new_parameter_points = p_box->SubDivide();

            // Evaluate at these points.
            RunEvaluationsForThesePoints(new_parameter_points, p_file);
        }
    }

    p_file->close();
}

template<unsigned DIM>
void LookupTableGenerator<DIM>::RunEvaluationsForThesePoints(CornerSet setOfPoints, out_stream& rFile)
{
    // Setup variables to control threading
    unsigned num_threads = setOfPoints.size();
    boost::scoped_array<ThreadInputData> thread_data(new ThreadInputData[num_threads]);
    //struct ThreadInputData thread_data[num_threads];
    boost::scoped_array<void*> answers(new void*[num_threads]);
    int return_code;

    // Generate the threads
    boost::scoped_array<pthread_t> threads(new pthread_t[num_threads]);

    // Create a couple of counters for convenience
    CornerSetIter iter;
    int i;

    /**
     * This loop launches each of the threads.
     */
    for (iter=setOfPoints.begin(), i=0;
         iter != setOfPoints.end();
         ++iter, ++i)
    {
        std::vector<double> scalings;
        for (unsigned j=0; j<DIM; j++)
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

        // Launch the ThreadedActionPotential method on this thread
        return_code = pthread_create(&threads[i],
                                     NULL,
                                     ThreadedActionPotential,
                                     (void*)&thread_data[i]);

        assert(0 == return_code); // Check launch worked OK
    }

    /*
     * This loop gets the answers back from all the threads.
     */
    for (iter = setOfPoints.begin(), i=0;
         iter != setOfPoints.end();
         ++iter, ++i)
    {
        // Get the answers back
        return_code = pthread_join(threads[i], &answers[i]);
        assert(0 == return_code);

        // Translate back from the structs to sensible formats.
        ThreadReturnData *thread_results = (ThreadReturnData *) answers[i];
        if (thread_results->exceptionOccurred)
        {
            EXCEPTION("A thread threw the exception: " << thread_results->exceptionMessage);
        }

        bool error_occurred = thread_results->errorOccurred;
        std::vector<double> results = thread_results->QoIs;
        c_vector<double, DIM>* p_scalings = *iter;
        delete thread_results; // Clean up memory

        // Store all the info in the master process and tell boxes about it.
        {
            boost::shared_ptr<ParameterPointData> data =
                    boost::shared_ptr<ParameterPointData>(new ParameterPointData(results, error_occurred));

            mParameterPoints.push_back(*p_scalings);
            mParameterPointData.push_back(data);
            mNumEvaluations++;

            // Tell all parameter boxes this information for future refinement.
            mpParentBox->AssignQoIValues(p_scalings, data);
            // This should have updated our error estimates in the ParameterPointData*

            std::stringstream line_of_output;
            line_of_output << std::setprecision(8);
            for (unsigned j=0; j<DIM; j++)
            {
                line_of_output << (*p_scalings)[j] << "\t";
            }
            line_of_output << error_occurred << "\t";
            for (unsigned j=0; j<results.size(); j++)
            {
                line_of_output << results[j] << "\t";
            }
            if (data->HasErrorEstimates())
            {
                unsigned num_estimates = data->rGetQoIErrorEstimates().size();
                line_of_output << num_estimates << "\t";
                for (unsigned j=0; j<num_estimates; j++)
                {
                    line_of_output << data->rGetQoIErrorEstimates()[j] << "\t";
                }
            }
            *rFile << line_of_output.str() << std::endl;
        }
    }
}



void* ThreadedActionPotential(void *argument)
{
   struct ThreadInputData *my_data;
   my_data = (struct ThreadInputData *) argument;

   std::vector<double> scalings = my_data->scalings;
   assert(scalings.size()==my_data->mParameterNames.size());

   SetupModel setup(my_data->mFrequency, my_data->mModelIndex); // Ten tusscher '06 at 1 Hz
   boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

   // Do parameter scalings
   for (unsigned i=0; i<scalings.size(); i++)
   {
       p_model->SetParameter(my_data->mParameterNames[i],
                             my_data->mUnscaledParameters[i]*(scalings[i]));
   }

   // Reset the state variables to the 'standard' steady state
   N_Vector state_vars = MakeNVector(my_data->mInitialConditions);
   p_model->SetStateVariables(state_vars);

   SingleActionPotentialPrediction ap_runner(p_model);
   ap_runner.SuppressOutput();
   ap_runner.SetMaxNumPaces(my_data->mMaxNumPaces);
   ap_runner.SetLackOfOneToOneCorrespondenceIsError();

   // Call the SingleActionPotentialPrediction methods.
   try
   {
       ap_runner.RunSteadyPacingExperiment();
   }
   catch (Exception& e)
   {
       ThreadReturnData* return_data = new ThreadReturnData;
       return_data->exceptionOccurred = true;
       return_data->exceptionMessage = e.GetShortMessage();

       DeleteVector(state_vars);
       pthread_exit(return_data);
   }

   bool error_occurred = ap_runner.DidErrorOccur();

   // Record the results
   std::vector<double> results;
   for (unsigned i=0; i<my_data->mQuantitiesToRecord.size(); i++)
   {
       if (error_occurred)
       {
           std::string error_code = ap_runner.GetErrorMessage();
           std::cout << "Lookup table generator reports that " << error_code << "\n" << std::flush;
           // We could use different numerical codes for different errors here if we wanted to.
           if (   (error_code=="NoActionPotential_2" || error_code =="NoActionPotential_3")
               && (my_data->mQuantitiesToRecord[i]==Apd90 || my_data->mQuantitiesToRecord[i]==Apd50) )
           {
               // For an APD calculation failure on repolarisation put in the stimulus period.
               double stim_period = boost::static_pointer_cast<RegularStimulus>(p_model->GetStimulusFunction())->GetPeriod();
               results.push_back(stim_period);
           }
           else
           {
               // For everything else (failure to depolarize "NoActionPotential_1") just put in zero for now.
               results.push_back(0.0);
           }
           continue;
       }

       double temp;
       if (my_data->mQuantitiesToRecord[i]==Apd90)
       {
           temp = ap_runner.GetApd90();
       }
       if (my_data->mQuantitiesToRecord[i]==Apd50)
       {
           temp = ap_runner.GetApd50();
       }
       if (my_data->mQuantitiesToRecord[i]==UpstrokeVelocity)
       {
           temp = ap_runner.GetUpstrokeVelocity();
       }
       if (my_data->mQuantitiesToRecord[i]==PeakVoltage)
       {
           temp = ap_runner.GetPeakVoltage();
       }
       results.push_back(temp);
   }

   ThreadReturnData* return_data = new ThreadReturnData;
   return_data->QoIs = results;
   return_data->errorOccurred = error_occurred;
   return_data->exceptionOccurred = false;

   DeleteVector(state_vars);
   pthread_exit(return_data);
}

template<unsigned DIM>
std::vector<c_vector<double, DIM> > LookupTableGenerator<DIM>::GetParameterPoints()
{
    return mParameterPoints;
}

template<unsigned DIM>
std::vector<std::vector<double> > LookupTableGenerator<DIM>::GetFunctionValues()
{
    std::vector<std::vector<double> > results;
    for (unsigned i=0; i<mParameterPointData.size(); i++)
    {
        results.push_back(mParameterPointData[i]->GetQoIs());
    }
    return results;
}

template<unsigned DIM>
void LookupTableGenerator<DIM>::SetParameterToScale(const std::string& rMetadataName,
                                                    const double& rMin,
                                                    const double& rMax)
{
    if (mParameterNames.size() == DIM)
    {
        EXCEPTION("All parameters have been defined already. You need to expand the dimension of your Lookup table generator.");
    }

    if (mGenerationHasBegun)
    {
        EXCEPTION("SetParameterToScale cannot be called after GenerateLookupTable.");
    }

    SetupModel setup(1.0, mModelIndex); // model at 1 Hz
    boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

    if (!p_model->HasParameter(rMetadataName))
    {
        // Not all of the models have a distinct fast I_to component.
        // In this case we look for the complete I_to current instead.
        if (rMetadataName=="membrane_fast_transient_outward_current_conductance" &&
            p_model->HasParameter("membrane_transient_outward_current_conductance") )
        {
            WARNING(p_model->GetSystemName() <<
                    " does not have 'membrane_fast_transient_outward_current_conductance' labelled, "
                    "using combined Ito (fast and slow) instead...");
            mParameterNames.push_back("membrane_transient_outward_current_conductance");
        }
        else
        {
            EXCEPTION(p_model->GetSystemName() << " does not have '" << rMetadataName <<
                      "' labelled, please tag it in the CellML file.");
        }
    }
    else
    {
        mParameterNames.push_back(rMetadataName);
    }
    mMinimums.push_back(rMin);
    mMaximums.push_back(rMax);
}

template<unsigned DIM>
void LookupTableGenerator<DIM>::AddQuantityOfInterest(QuantityOfInterest quantity, double tolerance)
{
    if (mGenerationHasBegun)
    {
        EXCEPTION("AddQuantityOfInterest cannot be called after GenerateLookupTable.");
    }
    mQuantitiesToRecord.push_back(quantity);
    mQoITolerances.push_back(tolerance);
}

template<unsigned DIM>
void LookupTableGenerator<DIM>::SetMaxNumEvaluations(const unsigned& rMaxNumEvals)
{
    mMaxNumEvaluations = rMaxNumEvals;
}

template<unsigned DIM>
void LookupTableGenerator<DIM>::SetMaxVariationInRefinement(const unsigned& rMaxRefinementDifference)
{
    mMaxRefinementDifference = rMaxRefinementDifference;
}

template<unsigned DIM>
std::vector<std::vector<double> > LookupTableGenerator<DIM>::Interpolate(const std::vector<c_vector<double, DIM> >& rParameterPoints)
{
    std::vector<std::vector<double> > interpolated_values;

    for (unsigned i=0; i<rParameterPoints.size(); i++)
    {
        interpolated_values.push_back(mpParentBox->InterpolateQoIsAt(rParameterPoints[i]));
    }

    return interpolated_values;
}

template<unsigned DIM>
unsigned LookupTableGenerator<DIM>::GetNumEvaluations()
{
    return mNumEvaluations;
}

template<unsigned DIM>
void LookupTableGenerator<DIM>::SetMaxNumPaces(unsigned numPaces)
{
    mMaxNumPaces = numPaces;
}

template<unsigned DIM>
unsigned LookupTableGenerator<DIM>::GetMaxNumPaces()
{
    return mMaxNumPaces;
}

template<unsigned DIM>
void LookupTableGenerator<DIM>::SetPacingFrequency(double frequency)
{
    mFrequency = frequency;
}

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LookupTableGenerator)

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////
template class LookupTableGenerator<1>;
template class LookupTableGenerator<2>;
template class LookupTableGenerator<3>;
template class LookupTableGenerator<4>;
template class LookupTableGenerator<5>;
// Just up to 5D for now, may need to be bigger eventually.

