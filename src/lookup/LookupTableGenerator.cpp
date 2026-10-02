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

#include <chrono>
#include <condition_variable>
#include <deque>
#include <iomanip> // for setprecision()
#include <map>
#include <mutex>
#include <thread>

#include "FileFinder.hpp"
#include "LookupTableGenerator.hpp"
#include "ParameterBox.hpp"
#include "SetupModel.hpp"
#include "SingleActionPotentialPrediction.hpp"

struct ThreadInputData
{
    std::vector<std::string> mParameterNames;
    std::vector<double> mUnscaledParameters;
    std::vector<QuantityOfInterest> mQuantitiesToRecord;
    std::vector<double> mInitialConditions;
    unsigned mMaxNumPaces;
    unsigned mModelIndex;
    double mFrequency;
    double mVoltageThreshold;
};

/**
 * Run an action potential simulation with some parameters scaled, and work out the QoIs.
 * Safe to call from several threads at once.
 *
 * @param rData  Details of the model and QoIs.
 * @param rScalings  The scaling factor to apply to each parameter.
 * @param rQoIs  Filled with the quantities of interest.
 * @param rErrorCode  Set to the error code from the AP evaluation (0 if no error).
 */
void EvaluateActionPotential(const ThreadInputData& rData,
                             const std::vector<double>& rScalings,
                             std::vector<double>& rQoIs,
                             unsigned& rErrorCode);

/**
 * Runs evaluations of QoIs, one std::thread per evaluation, and hands back
 * the results in the order that they finish.
 *
 * Only the thread that owns this object should call its methods. The destructor waits for any
 * evaluations that are still running.
 */
class EvaluationThreads
{
public:
    /** The result of one evaluation */
    struct Result
    {
        /** The ID the evaluation was launched with */
        unsigned id;
        /** The QoIs that were evaluated */
        std::vector<double> qois;
        /** The error code from the evaluation */
        unsigned errorCode;
        /** Whether the evaluation threw an exception */
        bool exceptionOccurred;
        /** The message of the exception, if there was one */
        std::string exceptionMessage;
    };

    /**
     * Constructor
     *
     * @param rEvaluate  The function to run on each thread.
     * @param launchDelay  A pause (in seconds) after launching each thread.
     */
    EvaluationThreads(const LookupTableEvaluationFunction& rEvaluate, double launchDelay)
            : mEvaluate(rEvaluate),
              mLaunchDelay(launchDelay)
    {
    }

    /** Destructor - waits for all running evaluations to finish */
    ~EvaluationThreads()
    {
        for (auto& r_thread : mThreads)
        {
            r_thread.second.join();
        }
    }

    /** @return The number of evaluations that are running or have finished but not been collected. */
    unsigned GetNumRunning() const
    {
        return mThreads.size();
    }

    /**
     * Start an evaluation on a new thread.
     *
     * @param id  An ID for this evaluation, which is returned with its result.
     * @param rScalings  The parameter scalings to evaluate at.
     */
    void Launch(unsigned id, const std::vector<double>& rScalings)
    {
        assert(mThreads.find(id) == mThreads.end());
        mThreads[id] = std::thread([this, id, rScalings]()
        {
            Result result;
            result.id = id;
            result.errorCode = 0u;
            result.exceptionOccurred = false;
            try
            {
                mEvaluate(rScalings, result.qois, result.errorCode);
            }
            catch (Exception& e)
            {
                result.exceptionOccurred = true;
                result.exceptionMessage = e.GetShortMessage();
            }
            catch (std::exception& e)
            {
                result.exceptionOccurred = true;
                result.exceptionMessage = e.what();
            }
            {
                std::lock_guard<std::mutex> lock(mMutex);
                mFinished.push_back(result);
            }
            mFinishedCondition.notify_one();
        });

        // Historically we have had seg. faults when launching AP simulation threads simultaneously,
        // so we stagger them a little.
        if (mLaunchDelay > 0.0)
        {
            std::this_thread::sleep_for(std::chrono::duration<double>(mLaunchDelay));
        }
    }

    /**
     * Wait for any running evaluation to finish.
     *
     * @return Its result.
     */
    Result WaitForResult()
    {
        assert(!mThreads.empty());
        Result result;
        {
            std::unique_lock<std::mutex> lock(mMutex);
            mFinishedCondition.wait(lock, [this]() { return !mFinished.empty(); });
            result = mFinished.front();
            mFinished.pop_front();
        }
        mThreads[result.id].join();
        mThreads.erase(result.id);
        return result;
    }

private:
    /** The function that evaluates QoIs */
    LookupTableEvaluationFunction mEvaluate;

    /** A pause (in seconds) after launching each thread */
    double mLaunchDelay;

    /** The threads that have been launched and not yet collected, by ID */
    std::map<unsigned, std::thread> mThreads;

    /** Protects #mFinished */
    std::mutex mMutex;

    /** Signalled when an evaluation finishes */
    std::condition_variable mFinishedCondition;

    /** Results of evaluations that have finished but not been collected */
    std::deque<Result> mFinished;
};

/**
 * @return The number of threads to use by default, the number of cores on this machine.
 */
unsigned DefaultNumThreads()
{
    unsigned num_cores = std::thread::hardware_concurrency();
    return (num_cores > 0u) ? num_cores : 1u;
}

/* Private constructor - just for archiving */
template <unsigned DIM>
LookupTableGenerator<DIM>::LookupTableGenerator()
    : AbstractUntemplatedLookupTableGenerator(),
      mModelIndex(0u),
      mpParentBox(NULL),
      mNumThreads(DefaultNumThreads()){};

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
      mVoltageThreshold(-50.0),
      mNumThreads(DefaultNumThreads())
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

    // Do a few special things the first time round (not needed for a test evaluation function).
    if (!mGenerationHasBegun && !mEvaluationFunctionForTesting)
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
    }

    // Work out how to evaluate QoIs at each point.
    LookupTableEvaluationFunction evaluate = mEvaluationFunctionForTesting;
    double launch_delay = 0.0;
    if (!evaluate)
    {
        ThreadInputData input_data;
        input_data.mParameterNames = mParameterNames;
        input_data.mUnscaledParameters = mUnscaledParameters;
        input_data.mQuantitiesToRecord = mQuantitiesToRecord;
        input_data.mInitialConditions = mInitialConditions;
        input_data.mMaxNumPaces = mMaxNumPaces;
        input_data.mModelIndex = mModelIndex;
        input_data.mFrequency = mFrequency;
        input_data.mVoltageThreshold = mVoltageThreshold;
        evaluate = [input_data](const std::vector<double>& rScalings, std::vector<double>& rQoIs, unsigned& rErrorCode)
        {
            EvaluateActionPotential(input_data, rScalings, rQoIs, rErrorCode);
        };
        launch_delay = 0.1; // seconds
    }

    if (!mGenerationHasBegun)
    {
        // Initial scalings
        CornerSet set_of_points = mpParentBox->GetCorners();
        assert(set_of_points.size() == pow(2, DIM));

        // Run these initial evaluations multi-threaded.
        RunEvaluationsForThesePoints(set_of_points, evaluate, launch_delay, p_file);

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
        bool meets_tolerance = RefineForQuantityOfInterest(quantitiy_idx, evaluate, launch_delay, p_file);

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
bool LookupTableGenerator<DIM>::QueuedBoxPriorityCompare::operator()(const QueuedBox& rA, const QueuedBox& rB) const
{
    if (rA.numErrorCodes != rB.numErrorCodes)
    {
        // We prioritise refining boxes in well-behaved space over those on edges of regions with errors.
        return rA.numErrorCodes < rB.numErrorCodes;
    }
    if (rA.errorEstimate != rB.errorEstimate)
    {
        return rA.errorEstimate > rB.errorEstimate;
    }
    if (rA.refinementLevel != rB.refinementLevel)
    {
        return rA.refinementLevel < rB.refinementLevel;
    }
    // Leaf boxes don't overlap, so no two have the same minimum corner.
    return c_vector_compare<DIM>()(&rA.min, &rB.min);
}

template <unsigned DIM>
bool LookupTableGenerator<DIM>::QueuedBoxLevelCompare::operator()(const QueuedBox& rA, const QueuedBox& rB) const
{
    if (rA.refinementLevel != rB.refinementLevel)
    {
        return rA.refinementLevel < rB.refinementLevel;
    }
    return QueuedBoxPriorityCompare()(rA, rB);
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::CollectLeafBoxes(ParameterBox<DIM>* pBox, std::vector<ParameterBox<DIM>*>& rLeaves)
{
    if (!pBox->mAmParent)
    {
        rLeaves.push_back(pBox);
        return;
    }
    for (unsigned i = 0; i < pBox->mDaughterBoxes.size(); i++)
    {
        CollectLeafBoxes(pBox->mDaughterBoxes[i], rLeaves);
    }
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::EnqueueIfNeedsRefinement(ParameterBox<DIM>* pBox,
                                                         unsigned quantityIndex,
                                                         std::set<QueuedBox, QueuedBoxPriorityCompare>& rQueue,
                                                         std::set<QueuedBox, QueuedBoxLevelCompare>& rQueueByLevel)
{
    if (!pBox->DoesBoxNeedFurtherRefinement(mQoITolerances[quantityIndex], quantityIndex))
    {
        return;
    }
    QueuedBox queued_box;
    queued_box.pBox = pBox;
    queued_box.numErrorCodes = pBox->GetNumErrors();
    queued_box.errorEstimate = pBox->GetMaxErrorInQoIEstimateInThisBox(quantityIndex);
    queued_box.refinementLevel = pBox->GetRefinementLevel();
    queued_box.min = pBox->mMin;
    rQueue.insert(queued_box);
    rQueueByLevel.insert(queued_box);
}

template <unsigned DIM>
bool LookupTableGenerator<DIM>::RefineForQuantityOfInterest(unsigned quantityIndex,
                                                            const LookupTableEvaluationFunction& rEvaluate,
                                                            double launchDelay,
                                                            out_stream& rFile)
{
    // A subdivision into 2^DIM boxes (the original refinement scheme) counts as DIM refinement levels,
    // so we scale the maximum difference in refinement to keep its meaning the same.
    unsigned max_refinement_level_difference = UNSIGNED_UNSET;
    if (mMaxRefinementDifference != UNSIGNED_UNSET)
    {
        max_refinement_level_difference = mMaxRefinementDifference * DIM;
    }

    // The boxes that need refining. All their corners have been evaluated, so their rankings won't change.
    std::set<QueuedBox, QueuedBoxPriorityCompare> queue;
    std::set<QueuedBox, QueuedBoxLevelCompare> queue_by_level;
    std::vector<ParameterBox<DIM>*> leaves;
    CollectLeafBoxes(mpParentBox, leaves);
    for (unsigned i = 0; i < leaves.size(); i++)
    {
        EnqueueIfNeedsRefinement(leaves[i], quantityIndex, queue, queue_by_level);
    }
    unsigned most_refined_level = mpParentBox->GetMostRefinedChild()->GetRefinementLevel();

    // Boxes that have been bisected, and are waiting for evaluations at the new points before their
    // daughters can be queued.
    struct Refinement
    {
        std::vector<ParameterBox<DIM>*> daughters;
        unsigned numPointsOutstanding;
    };
    std::vector<Refinement> refinements;

    // Points that are being evaluated (or waiting for a thread), and the refinements waiting for each.
    // A point can be shared by refinements of neighbouring boxes, but it is only evaluated once.
    std::map<c_vector<double, DIM>*, std::vector<unsigned>, c_vector_compare<DIM> > points_in_progress;
    std::deque<c_vector<double, DIM>*> points_to_launch;
    std::map<unsigned, c_vector<double, DIM>*> launched_points;
    unsigned next_launch_id = 0u;

    EvaluationThreads threads(rEvaluate, launchDelay);

    while (true)
    {
        // Refine boxes from the top of the queue until there is enough work for all the threads.
        while (threads.GetNumRunning() + points_to_launch.size() < mNumThreads
               && mNumEvaluations + points_in_progress.size() < mMaxNumEvaluations
               && !queue.empty())
        {
            QueuedBox next_box = *(queue.begin());

            // Check the selected box isn't going to refine one area too much,
            // if it is refine the least refined area instead.
            const QueuedBox& r_least_refined = *(queue_by_level.begin());
            if (max_refinement_level_difference != UNSIGNED_UNSET
                && most_refined_level - r_least_refined.refinementLevel >= max_refinement_level_difference
                && next_box.refinementLevel == most_refined_level)
            {
                next_box = r_least_refined;
            }
            queue.erase(next_box);
            queue_by_level.erase(next_box);

            // Bisect this box along the dimension with the largest error estimate
            // (NB if we GetCorners() after this, it includes the new points and makes no sense!).
            ParameterBox<DIM>* p_box = next_box.pBox;
            unsigned dimension_to_split = p_box->ChooseDimensionToSplit(quantityIndex);
            CornerSet new_points = p_box->SubDivide(dimension_to_split);
            most_refined_level = std::max(most_refined_level, next_box.refinementLevel + 1u);

            Refinement refinement;
            refinement.daughters = p_box->mDaughterBoxes;
            refinement.numPointsOutstanding = new_points.size();
            refinements.push_back(refinement);
            const unsigned refinement_idx = refinements.size() - 1u;

            for (CornerSetIter iter = new_points.begin(); iter != new_points.end(); ++iter)
            {
                if (points_in_progress.find(*iter) == points_in_progress.end())
                {
                    points_to_launch.push_back(*iter);
                }
                points_in_progress[*iter].push_back(refinement_idx);
            }

            // If the new points had all been evaluated already, the daughters have error estimates now.
            if (new_points.empty())
            {
                for (unsigned i = 0; i < refinement.daughters.size(); i++)
                {
                    EnqueueIfNeedsRefinement(refinement.daughters[i], quantityIndex, queue, queue_by_level);
                }
            }
        }

        // Start evaluating points if there are threads free.
        while (!points_to_launch.empty() && threads.GetNumRunning() < mNumThreads)
        {
            c_vector<double, DIM>* p_point = points_to_launch.front();
            points_to_launch.pop_front();
            std::vector<double> scalings(p_point->begin(), p_point->end());
            launched_points[next_launch_id] = p_point;
            threads.Launch(next_launch_id, scalings);
            next_launch_id++;
        }

        if (threads.GetNumRunning() == 0u)
        {
            // Nothing is running, so either the queue is empty or we have done enough evaluations.
            assert(points_to_launch.empty());
            break;
        }

        // Wait for an evaluation to finish and record it.
        EvaluationThreads::Result result = threads.WaitForResult();
        if (result.exceptionOccurred)
        {
            EXCEPTION("A thread threw the exception: " << result.exceptionMessage);
        }
        c_vector<double, DIM>* p_point = launched_points[result.id];
        launched_points.erase(result.id);
        RecordEvaluation(p_point, result.qois, result.errorCode, rFile);

        // Queue the daughters of any refinements that now have all their new points evaluated.
        const std::vector<unsigned>& r_waiting_refinements = points_in_progress[p_point];
        for (unsigned i = 0; i < r_waiting_refinements.size(); i++)
        {
            Refinement& r_refinement = refinements[r_waiting_refinements[i]];
            assert(r_refinement.numPointsOutstanding > 0u);
            r_refinement.numPointsOutstanding--;
            if (r_refinement.numPointsOutstanding == 0u)
            {
                for (unsigned j = 0; j < r_refinement.daughters.size(); j++)
                {
                    assert(r_refinement.daughters[j]->mAllCornersEvaluated);
                    EnqueueIfNeedsRefinement(r_refinement.daughters[j], quantityIndex, queue, queue_by_level);
                }
            }
        }
        points_in_progress.erase(p_point);
    }

    if (queue.empty())
    {
        std::cout << "Error estimates are within requested tolerances... finishing\n"
                  << std::flush;
        return true;
    }
    return false;
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::RecordEvaluation(c_vector<double, DIM>* pPoint,
                                                 const std::vector<double>& rQoIs,
                                                 unsigned errorCode,
                                                 out_stream& rFile)
{
    // Store all the info in the master process and tell boxes about it.
    boost::shared_ptr<ParameterPointData> data = boost::shared_ptr<ParameterPointData>(
        new ParameterPointData(rQoIs, errorCode));

    mParameterPoints.push_back(*pPoint);
    mParameterPointData.push_back(data);
    mNumEvaluations++;

    // Tell all parameter boxes this information for future refinement.
    mpParentBox->AssignQoIValues(pPoint, data);
    // This should have updated our error estimates in the ParameterPointData*

    std::stringstream line_of_output;
    line_of_output << std::setprecision(8);
    for (unsigned j = 0; j < DIM; j++)
    {
        line_of_output << (*pPoint)[j] << "\t";
    }
    line_of_output << errorCode;
    for (unsigned j = 0; j < rQoIs.size(); j++)
    {
        line_of_output << "\t" << rQoIs[j];
    }
    if (data->HasErrorEstimates())
    {
        unsigned num_estimates = data->rGetQoIErrorEstimates().size();
        line_of_output << "\t" << num_estimates;
        for (unsigned j = 0; j < num_estimates; j++)
        {
            line_of_output << "\t" << data->rGetQoIErrorEstimates()[j];
            // A report on progress towards meeting the tolerance on this QoI.
            line_of_output << "\t" << mpParentBox->ReportPercentageOfSpaceWhereToleranceIsMetForQoI(mQoITolerances[j], j);
        }
    }

    *rFile << line_of_output.str() << std::endl;
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::RunEvaluationsForThesePoints(
    CornerSet setOfPoints,
    const LookupTableEvaluationFunction& rEvaluate,
    double launchDelay,
    out_stream& rFile)
{
    std::vector<c_vector<double, DIM>*> points(setOfPoints.begin(), setOfPoints.end());
    std::vector<EvaluationThreads::Result> results(points.size());

    // Run up to mNumThreads evaluations at once.
    {
        EvaluationThreads threads(rEvaluate, launchDelay);
        unsigned num_launched = 0u;
        unsigned num_finished = 0u;
        while (num_finished < points.size())
        {
            while (num_launched < points.size() && threads.GetNumRunning() < mNumThreads)
            {
                std::vector<double> scalings(points[num_launched]->begin(), points[num_launched]->end());
                threads.Launch(num_launched, scalings);
                num_launched++;
            }
            EvaluationThreads::Result result = threads.WaitForResult();
            if (result.exceptionOccurred)
            {
                EXCEPTION("A thread threw the exception: " << result.exceptionMessage);
            }
            results[result.id] = result;
            num_finished++;
        }
    }

    // Record the results in the order of the set, so output doesn't depend on thread timings.
    for (unsigned i = 0; i < points.size(); i++)
    {
        RecordEvaluation(points[i], results[i].qois, results[i].errorCode, rFile);
    }
}

void EvaluateActionPotential(const ThreadInputData& rData,
                             const std::vector<double>& rScalings,
                             std::vector<double>& rQoIs,
                             unsigned& rErrorCode)
{
    assert(rScalings.size() == rData.mParameterNames.size());

    SetupModel setup(rData.mFrequency,
                     rData.mModelIndex); // Ten tusscher '06 at 1 Hz
    boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();

    // Do parameter scalings
    for (unsigned i = 0; i < rScalings.size(); i++)
    {
        std::string param_name;
        if (p_model->HasParameter(rData.mParameterNames[i]))
        {
            param_name = rData.mParameterNames[i];
        }
        else
        {
            param_name = rData.mParameterNames[i] + "_scaling_factor";
        }
        p_model->SetParameter(param_name,
                              rData.mUnscaledParameters[i] * (rScalings[i]));
    }

    // Reset the state variables to the 'standard' steady state
    p_model->SetStateVariables(rData.mInitialConditions);

    SingleActionPotentialPrediction ap_runner(p_model);
    ap_runner.SuppressOutput();
    ap_runner.SetMaxNumPaces(rData.mMaxNumPaces);
    ap_runner.SetLackOfOneToOneCorrespondenceIsError();
    ap_runner.SetVoltageThresholdForRecordingAsActionPotential(
        rData.mVoltageThreshold);

    // Call the SingleActionPotentialPrediction methods (any exception is passed back to the main thread).
    ap_runner.RunSteadyPacingExperiment();

    rErrorCode = ap_runner.GetErrorCode(); // 0 if there was no error

    // Record the results
    rQoIs.clear();
    for (unsigned i = 0; i < rData.mQuantitiesToRecord.size(); i++)
    {
        if (ap_runner.DidErrorOccur())
        {
            std::string error_message = ap_runner.GetErrorMessage();
            std::cout << "Lookup table generator reports that " << error_message
                      << "\n"
                      << std::flush;

            // We could use different numerical codes for different errors here if we
            // wanted to, but for QNet all AP errors are just set to -DBL_MAX.
            if (rData.mQuantitiesToRecord[i] == QNet)
            {
                rQoIs.push_back(-DBL_MAX);
                continue;
            }

            // We could use different numerical codes for different errors here if we
            // wanted to.
            if ((error_message == "NoActionPotential_2" || error_message == "NoActionPotential_3") && (rData.mQuantitiesToRecord[i] == Apd90 || rData.mQuantitiesToRecord[i] == Apd50))
            {
                // For an APD calculation failure on repolarisation put in the stimulus
                // period.
                double stim_period = boost::static_pointer_cast<RegularStimulus>(
                                         p_model->GetStimulusFunction())
                                         ->GetPeriod();
                rQoIs.push_back(stim_period);
            }
            else
            {
                // For everything else (failure to depolarize "NoActionPotential_1")
                // just put in zero for now.
                rQoIs.push_back(0.0);
            }
            continue;
        }

        // No error cases
        double temp;
        if (rData.mQuantitiesToRecord[i] == Apd90)
        {
            temp = ap_runner.GetApd90();
        }
        else if (rData.mQuantitiesToRecord[i] == Apd50)
        {
            temp = ap_runner.GetApd50();
        }
        else if (rData.mQuantitiesToRecord[i] == UpstrokeVelocity)
        {
            temp = ap_runner.GetUpstrokeVelocity();
        }
        else if (rData.mQuantitiesToRecord[i] == PeakVoltage)
        {
            temp = ap_runner.GetPeakVoltage();
        }
        else if (rData.mQuantitiesToRecord[i] == QNet)
        {
            temp = ap_runner.CalculateQNet();
        }
        rQoIs.push_back(temp);
    }
}

template <unsigned DIM>
void LookupTableGenerator<DIM>::SetNumThreads(unsigned numThreads)
{
    if (numThreads == 0u)
    {
        EXCEPTION("The number of threads must be at least one.");
    }
    mNumThreads = numThreads;
}

template <unsigned DIM>
unsigned LookupTableGenerator<DIM>::GetNumThreads() const
{
    return mNumThreads;
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
