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

#include <boost/assign.hpp>

#include "CheckpointArchiveTypes.hpp"

#include "Exception.hpp"
#include "Warnings.hpp"
#include "SetupModel.hpp"

#include "AbstractIvpOdeSolver.hpp"
#include "CellMLLoader.hpp"
#include "CommandLineArguments.hpp"
#include "FileFinder.hpp"
#include "RegularStimulus.hpp"

/* Mapping to implemnt backward compatibility with old hard coded model numbers */
const std::map<std::string, std::string> SetupModel::modelMapping = {{"1", "shannon_wang_puglisi_weber_bers_2004"},
                                                                     {"2", "ten_tusscher_model_2006_epi"},
                                                                     {"3", "mahajan_shiferaw_2008"},
                                                                     {"4", "hund_rudy_2004"},
                                                                     {"5", "grandi_pasqualini_bers_2010_ss"},
                                                                     {"6", "ohara_rudy_2011_endo"},
                                                                     {"7", "paci_hyttinen_aaltosetala_severi_ventricularVersion"},
                                                                     {"8", "ohara_rudy_cipa_v1_2017"}};

const std::unordered_set<std::string> SetupModel::forceNumericalJModels = {"hund_rudy_2004"};


SetupModel::SetupModel(const double& rHertz,
    unsigned modelIndex,
    boost::shared_ptr<OutputFileHandler> pHandler)
    : mpHandler(pHandler)
{
    /// Cvode cells use a CVODE solver regardless of which standard solver is passed in.
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    // Exceptions for wrong argument combinations
    if (modelIndex == UNSIGNED_UNSET && !CommandLineArguments::Instance()->OptionExists("--cellml") && !CommandLineArguments::Instance()->OptionExists("--model"))
    {
        EXCEPTION("Argument \"--model <index>\" is required");
    }
    if (modelIndex == UNSIGNED_UNSET && CommandLineArguments::Instance()->OptionExists("--cellml") && CommandLineArguments::Instance()->OptionExists("--model"))
    {
        EXCEPTION("You can only call ApPredict with the option '--model' OR '--cellml <file>'.");
    }

    // Figure out which cellml we need name
    std::string modelName;
    if (modelIndex != UNSIGNED_UNSET) //passed a number
    {
        modelName = std::to_string(modelIndex);
    }
    else if(CommandLineArguments::Instance()->OptionExists("--cellml"))
    {
        // passed a file name via --cellml
        WARNING("Argument --cellml <file> is deprecated use --model <file> instead.");
        modelName = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--cellml");
    }
    else
    {   
        // passed a an index, model name or file name via --cellml
        modelName = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--model");
    }

    // check if we have been given a file name and if so use that
    FileFinder cellml_file(modelName, RelativeTo::AbsoluteOrCwd);
    if (cellml_file.Exists())
    {
        if (mpHandler == NULL)
        {
            EXCEPTION("Trying to set up a dynamically loaded model without a working directory in SetupModel constructor.");
        }
        CellMLLoader loader(cellml_file, *mpHandler, {});
        mpModel = loader.LoadCvodeCell();
    }
    else
    { // we have been given a model name or number
        if(CommandLineArguments::Instance()->OptionExists("--cellml"))
        {
            EXCEPTION("Invalid file given with --cellml argument: " + modelName);
        }
        
        // Check if we have been given an index that can be mapped to a model name
        auto mapIterator = SetupModel::modelMapping.find(modelName);
        if(mapIterator != SetupModel::modelMapping.end())
        {
            modelName = mapIterator->second;
        }

        // Create model using factory
        if(ModelFactory::Exists(modelName , "AnalyticCvode"))
        {
            mpModel.reset((AbstractCvodeCell*)ModelFactory::Create(modelName , "AnalyticCvode", p_solver, p_stimulus));
        }
        else
        {  // throw an error if the model isn't found
            EXCEPTION("No model matches this index: " + modelName);
        }

        //set numerical Jacobean if needed
        mpModel->ForceUseOfNumericalJacobian(SetupModel::forceNumericalJModels.find(modelName) != SetupModel::forceNumericalJModels.end());
    }
    //std::cout << "* model = " << mpModel->GetSystemName() << std::endl;

    double s_magnitude = -15; // We will attempt to overwrite these with model specific ones below
    double s_duration = 3.0; // We will attempt to overwrite these with model specific ones below
    double s1_period = 1000.0 / rHertz; // ms - we may overwrite this with a model specific one if it is self exciting.

    // Use the default CellML stimulus amplitude and duration, but set start time and period to what we want.
    if (mpModel->HasCellMLDefaultStimulus())
    {
        boost::shared_ptr<RegularStimulus> p_reg_stim = mpModel->UseCellMLDefaultStimulus();
        s_magnitude = p_reg_stim->GetMagnitude();
        s_duration = p_reg_stim->GetDuration();
    }
    else if (mpModel->HasAttribute("SuggestedCycleLength"))
    {
        // If the model has no stimulus current it can be given this attribute tag
        // so we know roughly what cycle length to use on it...
        s1_period = mpModel->GetAttribute("SuggestedCycleLength");
        //std::cout << "s1 period = " << s1_period << std::endl;
    }

    if (CommandLineArguments::Instance()->OptionExists("--pacing-stim-duration"))
    {
        s_duration = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--pacing-stim-duration");
    }

    if (CommandLineArguments::Instance()->OptionExists("--pacing-stim-magnitude"))
    {
        s_magnitude = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--pacing-stim-magnitude");
    }

    // If this is the qNet case, load up some sensible steady state values:
    // Load archive of 0.5Hz steady state variables.
    if (modelIndex == 8u && fabs(s1_period - 2000.0) < 1e-4)
    {
        FileFinder archive_file("projects/ApPredict/test/data/ord_cipa_0.5Hz_state_vars.arch", RelativeTo::ChasteSourceRoot);
        if (archive_file.IsFile())
        {
            std::string archive_filename = archive_file.GetAbsolutePath();
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            std::vector<double> state_vars;
            input_arch >> state_vars;
            mpModel->SetStateVariables(state_vars);
        }
    }

    // We always use this so graphs look nice.
    double s_start = 1.0; // ms
    boost::shared_ptr<RegularStimulus> p_regular_stimulus(new RegularStimulus(s_magnitude, s_duration, s1_period, s_start));

    mpModel->SetStimulusFunction(p_regular_stimulus); // Assign the regular stimulus to the cell's stimulus
    mpModel->SetTolerances(1e-8, 1e-8);
}

boost::shared_ptr<AbstractCvodeCell> SetupModel::GetModel()
{
    return mpModel;
}
