/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "SetupModel.hpp"
#include "Exception.hpp"

#include "RegularStimulus.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "CommandLineArguments.hpp"
#include "CellMLLoader.hpp"
#include "FileFinder.hpp"

#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
//#include "Shannon2004Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "hund_rudy_2004Cvode.hpp"
#include "mahajan_shiferaw_2008Cvode.hpp"
#include "grandi_pasqualini_bers_2010_ssCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "paci_hyttinen_aaltosetala_severi_ventricularVersionCvode.hpp"

SetupModel::SetupModel(const double& rHertz,
                       unsigned modelIndex,
                       boost::shared_ptr<OutputFileHandler> pHandler)
 : mpHandler(pHandler)
{
    /// Cvode cells use a CVODE solver regardless of which standard solver is passed in.
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;


    // If modelIndex is specified, then we have to use that, and ignore command line.
    if (modelIndex == UNSIGNED_UNSET && CommandLineArguments::Instance()->OptionExists("--cellml"))
    {
        // Try to use a dynamically loaded model
        if(CommandLineArguments::Instance()->OptionExists("--model"))
        {
            EXCEPTION("You can only call ApPredict with the option '--model' OR '--cellml <file>'.");
        }
        std::string cellml_file_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--cellml");
        FileFinder cellml_file(cellml_file_path, RelativeTo::CWD);
        std::vector<std::string> options = boost::assign::list_of("--expose-annotated-variables");
        if (mpHandler==NULL)
        {
            EXCEPTION("Trying to set up a dynamically loaded model without a working directory in SetupModel constructor.");
        }
        CellMLLoader loader(cellml_file, *mpHandler, options);
        mpModel = loader.LoadCvodeCell();
    }
    else // Using a hardcoded model
    {
        if (modelIndex == UNSIGNED_UNSET)
        {
            if (!CommandLineArguments::Instance()->OptionExists("--model"))
            {
                EXCEPTION("Argument \"--model <index>\" is required");
            }
            modelIndex = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--model");
        }
        switch (modelIndex)
        {
            case 1u:
                // This one is from the cellml project - more metadata.
                mpModel.reset(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus));
                // This one is from the Chaste source
                //mpModel.reset(new CellShannon2004FromCellMLCvode(p_solver, p_stimulus));
                break;
            case 2u:
                mpModel.reset(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus));
                break;
            case 3u:
                mpModel.reset(new Cellmahajan_shiferaw_2008FromCellMLCvode(p_solver, p_stimulus));
                break;
            case 4u:
                mpModel.reset(new Cellhund_rudy_2004FromCellMLCvode(p_solver, p_stimulus));
                break;
            case 5u:
                mpModel.reset(new Cellgrandi_pasqualini_bers_2010_ssFromCellMLCvode(p_solver, p_stimulus));
                break;
            case 6u:
                mpModel.reset(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus));
                break;
            case 7u:
                mpModel.reset(new Cellpaci_hyttinen_aaltosetala_severi_ventricularVersionFromCellMLCvode(p_solver, p_stimulus));
                break;
            default:
                EXCEPTION("No model matches this index");
        }
    }
    //std::cout << "* model = " << mpModel->GetSystemName() << "\n";

    double s_magnitude = -15; // We will attempt to overwrite these with model specific ones below
    double s_duration = 3.0; // We will attempt to overwrite these with model specific ones below
    double s1_period = 1000.0/rHertz;  // ms - we may overwrite this with a model specific one if it is self exciting.

    // Use the default CellML stimulus amplitude and duration, but set start time and period to what we want.
    if (mpModel->HasCellMLDefaultStimulus())
    {
        boost::shared_ptr<RegularStimulus> p_reg_stim = mpModel->UseCellMLDefaultStimulus();
        s_magnitude = p_reg_stim->GetMagnitude();
        s_duration= p_reg_stim->GetDuration();
    }
    else if(mpModel->HasAttribute("SuggestedCycleLength"))
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


    // We always use this so graphs look nice.
    double s_start = 1.0;            // ms
    boost::shared_ptr<RegularStimulus> p_regular_stimulus(new RegularStimulus(s_magnitude, s_duration, s1_period, s_start));

    mpModel->SetStimulusFunction(p_regular_stimulus); // Assign the regular stimulus to the cell's stimulus
    mpModel->SetTolerances(1e-6, 1e-8);
}

boost::shared_ptr<AbstractCvodeCell> SetupModel::GetModel()
{
    return mpModel;
}
