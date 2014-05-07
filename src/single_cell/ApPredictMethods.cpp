/*

Copyright (c) 2005-2014, University of Oxford.
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

#include <fstream>
#include <numeric> // for std::accumulate
#include <sys/stat.h> // For system commands to download and unpack Lookup Table file.

// Chaste source includes
#include "CheckpointArchiveTypes.hpp"

#include "ApPredictMethods.hpp"

#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "CommandLineArguments.hpp"
#include "FileFinder.hpp"

#include "AbstractDataStructure.hpp"
#include "DoseCalculator.hpp"
#include "ActionPotentialDownsampler.hpp"
#include "BayesianInferer.hpp"

#include "ZeroStimulus.hpp"
#include "RegularStimulus.hpp"

#include "SetupModel.hpp"
#include "SteadyStateRunner.hpp"
#include "ProgressReporter.hpp"

std::string ApPredictMethods::PrintArguments()
{
	std::string message = "ApPredict::Please provide these inputs:\n" + SetupModel::PrintArguments();
	message += PrintCommonArguments();
	return message;
}

std::string ApPredictMethods::PrintCommonArguments()
{
    std::string message = "*\n"
            "* Specifying pacing frequency:\n"
            "* --pacing-freq       Pacing frequency (Hz) (optional - defaults to 1Hz)\n"
            "* --pacing-max-time   Maximum time for which to pace the cell model in MINUTES\n"
            "*                     (optional - defaults to time for 10,000 paces at this frequency)\n" // Set in AbstractSteadyStateRunner constructor!
            "*\n"
            "* SPECIFYING DRUG PROPERTIES dose-response properties for each channel:\n"
            "* Channels are named:\n"
            "* * herg (IKr current - hERG),\n"
            "* * na (fast sodium current - NaV1.5),\n"
            "* * cal (L-type calcium current- CaV1.2),\n"
            "* * iks (IKs current - KCNQ1 + MinK),\n"
            "* * ik1 (IK1 current - KCNN4 a.k.a. KCa3.1),\n"
            "* * ito ([fast] Ito current - Kv4.3 + KChIP2.2).\n"
            "*\n"
            "For each channel you specify dose-response parameters [multiple entries for repeat experiments]\n"
            "*   EITHER with IC50 values (in uM), for example for 'hERG':\n"
            "* --ic50-herg     hERG IC50    (optional - defaults to \"no affect\")\n"
            "*   OR with pIC50 values (in log M):\n"
            "* --pic50-herg    hERG pIC50   (optional - defaults to \"no affect\")\n"
            "*     (you can use a mixture of these for different channels if you wish, \n"
            "*     e.g. --ic50-herg 16600 --pic50-na 5.3 )\n"
            "*   AND specify Hill coefficients (dimensionless):\n"
            "* --hill-herg     hERG Hill    (optional - defaults to \"1.0\")\n"
            "*\n"
            "* SPECIFYING CONCENTRATIONS:\n"
            "* --plasma-concs  A list of plasma concentrations at which to test (uM)\n"
            "* OR alternatively:\n"
            "* --plasma-conc-high  Highest plasma concentration to test (uM)\n"
            "* --plasma-conc-low   Lowest  plasma concentration to test (uM) \n"
            "*                                  (optional - defaults to 0)\n"
            "*\n"
            "* both ways of specifying test concentrations have the following optional arguments\n"
            "* --plasma-conc-count  Number of intermediate plasma concentrations to test \n"
            "*                                  (optional - defaults to 0 (for --plasma-concs) or 11 (for --plasma-conc-high)\n"
            "* --plasma-conc-logscale  Use a log spacing for the plasma concentrations \n"
            "*                                  (optional - defaults to false)\n"
            "*\n"
            "* UNCERTAINTY QUANTIFICATION:\n"
            "* --credible-intervals  This flag must be present.\n"
            "* Then to specify 'spread' parameters for assay variability - for use with Lookup Tables:\n"
            "* --pic50-spread-herg      (for each channel that you are providing ic50/pic50 values for,\n"
            "* --hill-spread-herg        herg is just given as an example)\n"
            "*   (for details of what these spread parameters are see 'sigma' and '1/beta' in Table 1 of:\n"
            "*    Elkins et al. 2013  Journal of Pharmacological and Toxicological \n"
            "*    Methods, 68(1), 112-122. doi: 10.1016/j.vascn.2013.04.007 )\n"
            "*\n";

    return message;
}

void ApPredictMethods::ReadInIC50AndHill(std::vector<double>& rIc50s,
                                         std::vector<double>& rHills,
                                         const unsigned channelIdx)
{
    const std::string channel = mShortNames[channelIdx];
	CommandLineArguments* p_args = CommandLineArguments::Instance();
	bool read_ic50s = false;
	bool read_hills = false;

	// Try loading any arguments given as IC50s
    if (p_args->OptionExists("--ic50-" + channel))
    {
        rIc50s = p_args->GetDoublesCorrespondingToOption("--ic50-" + channel);
        if(p_args->OptionExists("--pic50-" + channel))
        {
        	EXCEPTION("Duplicate arguments, you cannot specify both IC50 and pIC50 for " << channel << " channel.");
        }
        read_ic50s = true;
    }
    // If those don't exist try loading pIC50s.
    else if(p_args->OptionExists("--pic50-" + channel))
    {
        rIc50s.clear();
        std::vector<double> pIC50s = p_args->GetDoublesCorrespondingToOption("--pic50-" + channel);
        for (unsigned i=0; i<pIC50s.size(); i++)
        {
            rIc50s.push_back(AbstractDataStructure::ConvertPic50ToIc50(pIC50s[i]));
        }
        read_ic50s = true;
    }

    // Try loading any Hills
    if (p_args->OptionExists("--hill-" + channel))
    {
        rHills = p_args->GetDoublesCorrespondingToOption("--hill-" + channel);
        // But these must correspond to IC50s.
        if (!(rHills.size()==rIc50s.size()))
        {
            EXCEPTION("If you enter Hill coefficients, there must be one corresponding to each [p]IC50 measurement.");
        }
        read_hills = true;
    }

    // Collect any spread parameter information that has been inputted.
    if(p_args->OptionExists("--pic50-spread-" + channel))
    {
        mPic50Spreads[channelIdx] = p_args->GetDoubleCorrespondingToOption("--pic50-spread-" + channel);
    }
    if(p_args->OptionExists("--hill-spread-" + channel))
    {
        mHillSpreads[channelIdx] = p_args->GetDoubleCorrespondingToOption("--hill-spread-" + channel);
    }

    if (mSuppressOutput)
    {
        return;
    }

    // Print the data we read in to screen for provenance tracking.
    std::cout << "* " << channel;
    if (read_ic50s)
    {
        std::cout << " IC50s = ";
        for (unsigned i=0; i<rIc50s.size(); i++)
        {
            std::cout << rIc50s[i] << " ";
        }
        std::cout << " uM, ";
        if (read_hills)
        {
            std::cout << "Hills = ";
            for (unsigned i=0; i<rHills.size(); i++)
            {
                std::cout << rHills[i] << " ";
            }
            std::cout << "\n";
        }
        else
        {
            std::cout << "Hills = 1.0 (default) \n";
        }
    }
    else
    {
        std::cout << ": no drug affect\n";
    }
}

void ApPredictMethods::ApplyDrugBlock(boost::shared_ptr<AbstractCvodeCell> pModel,
                                      unsigned channel_index,
                                      const double default_conductance,
                                      const double concentration,
                                      const double iC50,
                                      const double hill)
{
    // Here we calculate the proportion of the different channels which are still active
    // (at this concentration of this drug)
    const double conductance_factor = AbstractDataStructure::CalculateConductanceFactor(concentration, iC50,  hill);

    // Some screen output for info.
    if (!mSuppressOutput) std::cout << "g_" << mShortNames[channel_index] <<  " factor = " << conductance_factor << "\n";// << std::flush;

    // Check the model has this parameter before we try and set it.
    if (pModel->HasParameter(mMetadataNames[channel_index]))
    {
        pModel->SetParameter(mMetadataNames[channel_index], default_conductance*conductance_factor);
    }
    else // We haven't got that conductance parameter, or at least it isn't labelled.
    {
        // If we aren't trying to change it - don't worry, just carry on.
        if (conductance_factor < 1)
        {
            // If the model hasn't got this channel conductance labelled,
            // (but we are trying to change it) throw an error.
            EXCEPTION(pModel->GetSystemName() << " does not have the current \"" << mMetadataNames[channel_index] << "\" labelled, but you have requested a block on this channel.");
        }
    }
}

ApPredictMethods::ApPredictMethods()
   : AbstractActionPotentialMethod(),
     mLookupTableAvailable(false),
     mComplete(false)
{
    // Here we list the possible drug blocks that can be applied with ApPredict
    mMetadataNames.push_back("membrane_fast_sodium_current_conductance");
    mShortNames.push_back("na");

    mMetadataNames.push_back("membrane_L_type_calcium_current_conductance");
    mShortNames.push_back("cal");

    mMetadataNames.push_back("membrane_rapid_delayed_rectifier_potassium_current_conductance");
    mShortNames.push_back("herg");

    mMetadataNames.push_back("membrane_slow_delayed_rectifier_potassium_current_conductance");
    mShortNames.push_back("iks");

    mMetadataNames.push_back("membrane_inward_rectifier_potassium_current_conductance");
    mShortNames.push_back("ik1");

    mMetadataNames.push_back("membrane_fast_transient_outward_current_conductance");
    mShortNames.push_back("ito");

    // We'll use DOUBLE_UNSET to start with for these spread parameters.
    for (unsigned i=0; i<6u; i++)
    {
        mPic50Spreads.push_back(DOUBLE_UNSET);
        mHillSpreads.push_back(DOUBLE_UNSET);
    }

    // There must be a 1:1 mapping between these...
    assert(mMetadataNames.size()==mShortNames.size());
}

void ApPredictMethods::SetUpLookupTables()
{
    CommandLineArguments* p_args = CommandLineArguments::Instance();

    if (!p_args->OptionExists("--credible-intervals"))
    {
        // The flag mLookupTableAvailable remains false, and we carry on as normal.
        return;
    }

    // We've only generated (up to) 4D lookup tables for now,
    // so don't bother if we are asking for ion channel blocks of other things
    if (p_args->OptionExists("--ic50-ik1")  ||
        p_args->OptionExists("--pic50-ik1") ||
        p_args->OptionExists("--ic50-ito")  ||
        p_args->OptionExists("--pic50-ito") )
    {
        EXCEPTION("Lookup table (for --credible-intervals) is currently only including IKr, IKs, INa and ICaL block, you have specified additional ones so quitting.");
    }

    // Here we will attempt to use any lookup table associated with this model and pacing rate.
    std::stringstream lookup_table_archive_name;
    lookup_table_archive_name << mpModel->GetSystemName() << "_4d_hERG_IKs_INa_ICaL_" <<  this->mHertz << "Hz_generator";

    // First see if there is a table available already in absolute or current working directory.
    FileFinder ascii_archive_file(lookup_table_archive_name.str() + ".arch", RelativeTo::AbsoluteOrCwd);
    FileFinder binary_archive_file(lookup_table_archive_name.str() + "_BINARY.arch", RelativeTo::AbsoluteOrCwd);

    // First we try loading the binary version of the archive, if it exists.
    if (binary_archive_file.IsFile())
    {
        std::cout << "Loading lookup table from binary archive into memory, this can take a few seconds..." << std::flush;
        double start = MPI_Wtime();

        // Create a pointer to the input archive
        std::ifstream ifs((binary_archive_file.GetAbsolutePath()).c_str(), std::ios::binary);
        boost::archive::binary_iarchive input_arch(ifs);

        // restore from the archive
        LookupTableGenerator<TABLE_DIM>* p_generator;
        input_arch >> p_generator;

        mpLookupTable.reset(p_generator);
        mLookupTableAvailable = true;

        double load_time = MPI_Wtime() - start;
        std::cout << " loaded in " << load_time << " secs.\nLookup table is available for generation of credible intervals.\n";

        // Since loading the binary archive works, we can try and get rid of the ascii one to clean up.
        if (ascii_archive_file.IsFile())
        {
            try
            {
                // The ascii file is not in a testoutput folder so we need to over-ride our usual safety checks.
                ascii_archive_file.DangerousRemove();
                std::cout << "Ascii lookup table archive file removed to tidy up, will use the binary one in future." << std::endl;
            }
            catch (Exception &e)
            {
                WARNING("Could not remove ascii lookup table archive, error was: " << e.GetMessage() << "\nSimulations continued anyway.");
            }
        }

        // We have finished and loaded the generator.
        return;
    }

    if (!ascii_archive_file.IsFile())
    {
        // If no archive exists, try to download and unpack one.
        std::string lookup_table_URL = "http://www.cs.ox.ac.uk/people/gary.mirams/files/" + lookup_table_archive_name.str() + ".arch.tgz";
        try
        {
            std::cout << "\n\nAttempting to download an action potential lookup table from:\n" << lookup_table_URL << "\n\n";
            EXPECT0(system, "wget " + lookup_table_URL);
            std::cout << "Download succeeded, unpacking...\n";
            EXPECT0(system, "tar xzf " + lookup_table_archive_name.str() + ".arch.tgz");
            std::cout << "Unpacking succeeded, removing .tgz file...\n";
            EXPECT0(system, "rm -f " + lookup_table_archive_name.str() + ".arch.tgz");
        }
        catch (Exception &e)
        {
            std::cout << "Could not download and unpack the Lookup Table archive, continuing without it...\n";
            return;
        }
    }

    // If there is no binary archive, then try loading the ascii version and creating a binary one for next time.
    if (ascii_archive_file.IsFile())
    {
        std::cout << "Loading lookup table from file into memory, this can take a few seconds..." << std::flush;
        double start = MPI_Wtime();

        // Create a pointer to the input archive
        std::ifstream ifs((ascii_archive_file.GetAbsolutePath()).c_str(), std::ios::binary);
        boost::archive::text_iarchive input_arch(ifs);

        // restore from the archive
        LookupTableGenerator<TABLE_DIM>* p_generator;
        input_arch >> p_generator;

        mpLookupTable.reset(p_generator);
        mLookupTableAvailable = true;

        double load_time = MPI_Wtime() - start;
        std::cout << " loaded in " << load_time << " secs.\nLookup table is available for generation of credible intervals.\n";

        try
        {
            std::cout << "Saving a binary version of the archive for faster loading next time..." << std::flush;
            // Save a binary version to speed things up next time round.
            LookupTableGenerator<TABLE_DIM>* const p_arch_generator = p_generator;
            std::ofstream binary_ofs(binary_archive_file.GetAbsolutePath().c_str(), std::ios::binary);
            boost::archive::binary_oarchive output_arch(binary_ofs);
            output_arch << p_arch_generator;
            std::cout << "done!\n";
        }
        catch (Exception &e)
        {
            WARNING("Did not manage to create binary lookup table archive. Error was: " << e.GetMessage() << "\nContinuing to use ascii archive.");
        }
    }
}

void ApPredictMethods::CalculateDoseResponseParameterSamples(const std::vector<std::vector<double> >& rIC50s,
                                                             const std::vector<std::vector<double> >& rHills)
{
    if (mLookupTableAvailable)
    {
        /*
         * Prepare an inferred set of IC50s and Hill coefficients
         * for use with the Lookup Table and credible interval calculations.
         */
        mSampledIc50s.resize(mMetadataNames.size());
        mSampledHills.resize(mMetadataNames.size());

        const unsigned num_samples = 1000u;

        // Work out vectors of inferred IC50 and Hills
        // Apply drug block on each channel
        for (unsigned channel_idx = 0 ; channel_idx<mMetadataNames.size(); channel_idx++)
        {
            // First just decide whether there is 'no effect here'.
            assert(rIC50s[channel_idx].size()>=1u);
            if (rIC50s[channel_idx].size()==1 && rIC50s[channel_idx][0]==-1)
            {
                // No effect here, so just map that out to all random runs
                for (unsigned i=0; i<num_samples; i++)
                {
                    mSampledIc50s[channel_idx].push_back(-1.0);
                    mSampledHills[channel_idx].push_back(-1.0);
                }
                // To next channel
                continue;
            }

            // Convert data to pIC50s
            std::vector<double> pIC50s;
            for (unsigned i=0; i<rIC50s[channel_idx].size(); i++)
            {
                pIC50s.push_back(AbstractDataStructure::ConvertIc50ToPic50(rIC50s[channel_idx][i]));
            }

            // Retrieve the Pic50 spread from the command line arguments.
            if (mPic50Spreads[channel_idx] == DOUBLE_UNSET)
            {
                EXCEPTION("No argument --pic50-spread-" << mShortNames[channel_idx] << " has been provided. Cannot calculate credible intervals without this.");
            }

            // Infer pIC50 spread.
            BayesianInferer ic50_inferer(PIC50);
            ic50_inferer.SetObservedData(pIC50s);
            ic50_inferer.SetSpreadOfUnderlyingDistribution(mPic50Spreads[channel_idx]);
            ic50_inferer.PerformInference();

            std::vector<double> inferred_pic50s = ic50_inferer.GetSampleMedianValue(num_samples); // Get 1000 inferred pIC50s
            for (unsigned i=0; i<num_samples; i++)
            {
                // Convert pIC50 back to IC50 and store it.
                mSampledIc50s[channel_idx].push_back(AbstractDataStructure::ConvertPic50ToIc50(inferred_pic50s[i]));
            }

            // If all Hill coefficient entries are positive, then we will use those for samples too.
            // This long-winded bit of code is just counting how many are positive (must be a better way!).
            bool all_hills_positive = false;
            unsigned temp_counter = 0u;
            for (unsigned i=0; i<rHills[channel_idx].size(); i++)
            {
                if (rHills[channel_idx][i] > 0.0)
                {
                    temp_counter++;
                }
            }
            if (temp_counter == rHills[channel_idx].size())
            {
                all_hills_positive = true;
            }

            if (all_hills_positive)
            {
                // Retrieve Hill spread parameters from stored Command line args.
                if (mHillSpreads[channel_idx] == DOUBLE_UNSET)
                {
                    EXCEPTION("No argument --hill-spread-" << mShortNames[channel_idx] << " has been provided. Cannot calculate credible intervals without this.");
                }

                // Infer Hill spread.
                BayesianInferer hill_inferer(HILL);
                hill_inferer.SetObservedData(rHills[channel_idx]);
                // This works with the Beta parameter, not the 1/Beta. So do 1/1/Beta to get Beta back!
                hill_inferer.SetSpreadOfUnderlyingDistribution(1.0/mHillSpreads[channel_idx]);
                hill_inferer.PerformInference();

                mSampledHills[channel_idx] = hill_inferer.GetSampleMedianValue(num_samples); // Get 1000 inferred Hills
            }
            else
            {
                // There have been not enough 'real' Hill coefficients specified,
                // so assume they are all missing (otherwise inference would go mad
                // with some positive and some negative).
                for (unsigned i=0; i<num_samples; i++)
                {
                    // We aren't going to attempt to do inference on Hills,
                    // just push back 'no measurement' for now.
                    mSampledHills[channel_idx].push_back(-1.0);
                }
            }
        }
    }
}

void ApPredictMethods::InterpolateFromLookupTableForThisConcentration(const unsigned conc_index)
{
    // If we don't have a lookup table, we aren't going to do confidence intervals.
    if (!mLookupTableAvailable)
    {
        return;
    }

    std::pair<double, double> credible_interval;

    // If this is the first concentration (control) say the percent change must be zero
    // or otherwise a small interpolation error will result
    // (from potentially running to different steady state with --pacing-max-time <x> ).
    if (conc_index==0u)
    {
        credible_interval.first = mApd90s[conc_index];
        credible_interval.second = mApd90s[conc_index];
        mApd90CredibleRegions[conc_index] = credible_interval;
        return;
    }

    // The first channel entry in mSampledIc50s will give us the number of random samples.
    const unsigned num_samples = mSampledIc50s[0].size();

    std::cout << "Calculating confidence intervals from Lookup Table...";
    std::vector<c_vector<double,4u> > sampling_points;
    for (unsigned rand_idx=0; rand_idx<num_samples; rand_idx++)
    {
        // In the lookup table the order of parameters is given in the filename:
        // "4d_hERG_IKs_INa_ICaL_generator.arch"
        c_vector<double,TABLE_DIM> sample_required_at;

        // IKr
        sample_required_at[0] = AbstractDataStructure::CalculateConductanceFactor(mConcs[conc_index],
                                                                                  mSampledIc50s[2][rand_idx],
                                                                                  mSampledHills[2][rand_idx]);
        // IKs
        sample_required_at[1] = AbstractDataStructure::CalculateConductanceFactor(mConcs[conc_index],
                                                                                  mSampledIc50s[3][rand_idx],
                                                                                  mSampledHills[3][rand_idx]);
        // INa
        sample_required_at[2] = AbstractDataStructure::CalculateConductanceFactor(mConcs[conc_index],
                                                                                  mSampledIc50s[0][rand_idx],
                                                                                  mSampledHills[0][rand_idx]);
        // ICaL
        sample_required_at[3] = AbstractDataStructure::CalculateConductanceFactor(mConcs[conc_index],
                                                                                  mSampledIc50s[1][rand_idx],
                                                                                  mSampledHills[1][rand_idx]);
        sampling_points.push_back(sample_required_at);
    }

    std::vector<std::vector<double> > predictions = mpLookupTable->Interpolate(sampling_points);
    assert(predictions.size()==mSampledIc50s[0].size());

    /*
     * Compile all the lookup table predictions into a vector that we can sort.
     */
    //std::vector<double> apd_50_predictions;
    std::vector<double> apd_90_predictions;
    for (unsigned rand_idx=0; rand_idx<num_samples; rand_idx++)
    {
        //apd_50_predictions.push_back(predictions[rand_idx][1]);
        apd_90_predictions.push_back(predictions[rand_idx][0]);
    }
    std::sort(apd_90_predictions.begin(), apd_90_predictions.end());

    // Now work out the confidence intervals, err on conservative side.
    double credible_region = 95.0; // %
    double tails = 0.01*(100 - credible_region)/2.0;
    credible_interval.first = apd_90_predictions[floor(tails*(double)(num_samples))];
    credible_interval.second = apd_90_predictions[ceil((1-tails)*(double)(num_samples))];

    mApd90CredibleRegions[conc_index] = credible_interval;
    std::cout << "done.\n";
}

void ApPredictMethods::Run()
{
    mProgramName = "Action Potential PreDiCT";
    mOutputFolder = "ApPredict_output/";

    SetupModel setup(this->mHertz); // This class will get model definition from command line.
    mpModel = setup.GetModel();

    SetUpLookupTables();

    CommonRunMethod();
}

void ApPredictMethods::CommonRunMethod()
{
    // Arguments that take default values
    std::vector<std::vector<double> > IC50s;
    std::vector<std::vector<double> > hills;

    std::vector<double> unset;
    unset.push_back(-1); // -1 is our code for 'unset'.

    // Set the initial values of these
    for (unsigned i=0; i<mMetadataNames.size(); i++)
    {
        IC50s.push_back(unset);
        hills.push_back(unset);
    }

    // Dose calculator asks for some arguments to do with plasma concentrations.
    DoseCalculator dose_calculator;
    mConcs = dose_calculator.GetConcentrations();

    // We check the desired parameters are present in the model, warn if not.
    // This method also changes some metadata names if the model has variants that will do, but aren't ideal
    // and warns if it does this.
    ParameterWrapper(mpModel, mMetadataNames);

    // Use a helper method to read in IC50 from either --ic50 or --pic50 arguments.
    // Note this is now in micro Molar (1x10^-6 Molar) as per most Pharma use.
    for (unsigned channel_idx=0; channel_idx<mMetadataNames.size(); channel_idx++)
    {
        ReadInIC50AndHill(IC50s[channel_idx], hills[channel_idx], channel_idx);
    }

    if (!mSuppressOutput)
    {
        std::cout << "* max free plasma concentration = " << mConcs.back()<< " uM\n"
                     "* min free plasma concentration = " << mConcs[0] << " uM\n"
                     "* number of plasma concentrations = " << mConcs.size() << "\n";// << std::flush;
    }

    // The following names are fixed and correspond to metadata tags.
    // We record the default parameter values that the model uses.
    // All the drug block models should include these parameter labels
    // But you only get a warning if not, so check the warnings...
    std::vector<double> default_conductances;
    for (unsigned channel_idx=0; channel_idx<mMetadataNames.size(); channel_idx++)
    {
        double default_value = 1.0;
        if (mpModel->HasParameter(mMetadataNames[channel_idx]))
        {
            default_value = mpModel->GetParameter(mMetadataNames[channel_idx]);
        }
        default_conductances.push_back(default_value);
    }

    CalculateDoseResponseParameterSamples(IC50s, hills);

    boost::shared_ptr<const AbstractOdeSystemInformation> p_ode_info = mpModel->GetSystemInformation();
    std::string model_name = mpModel->GetSystemName();
    boost::static_pointer_cast<RegularStimulus>(mpModel->GetStimulusFunction())->SetStartTime(5.0);

    // Set up foldernames for each model and protocol set.
    std::string foldername = mOutputFolder;

    // Make and clean the above directories.
    mpFileHandler.reset(new OutputFileHandler(foldername));

    // Print out a progress file for monitoring purposes.
    ProgressReporter progress_reporter(foldername, 0.0, (double)(mConcs.size()));
    progress_reporter.PrintInitialising();

    // Open files and write headers
    out_stream steady_voltage_results_file_html = mpFileHandler->OpenOutputFile("voltage_results.html");

    out_stream steady_voltage_results_file = mpFileHandler->OpenOutputFile("voltage_results.dat");
    *steady_voltage_results_file << "Concentration(uM)\tUpstrokeVelocity(mV/ms)\tPeakVm(mV)\tAPD50(ms)\tAPD90(ms)\t";
    if (mLookupTableAvailable)
    {
        *steady_voltage_results_file << "delta_APD90_lower(%),delta_APD90(%),delta_APD90_upper(%)\n";
    }
    else
    {
        *steady_voltage_results_file << "delta_APD90(%)\n";
    }

    *steady_voltage_results_file_html << "<html>\n<head><title>" << mProgramName << " results</title></head>\n";
    *steady_voltage_results_file_html << "<STYLE TYPE=\"text/css\">\n<!--\nTD{font-size: 12px;}\n--->\n</STYLE>\n";
    *steady_voltage_results_file_html << "<body>\n";
    *steady_voltage_results_file_html << "<table width=\"60%\" style=\"background-color:white\" border=\"1\" cellpadding=\"2\" cellspacing=\"0\">\n";
    *steady_voltage_results_file_html << "<tr><td>Concentration (uM)</td><td>Upstroke Velocity (mV/ms)</td><td>Peak membrane voltage (mV)</td><td>APD50 (ms)</td><td>APD90 (ms)</td><td>Change in APD90 (%)</td></tr>\n"; // Header line

    /**
     * START LOOP OVER EACH CONCENTRATION TO TEST WITH
     */
    mApd90CredibleRegions.resize(mConcs.size());
    double control_apd90 = 0;
    for (unsigned conc_index=0; conc_index<mConcs.size(); conc_index++)
    {
        progress_reporter.Update((double)(conc_index));
        std::cout << "Drug Conc = " << mConcs[conc_index] << " uM\n" ;//<< std::flush;

        // Apply drug block on each channel
        for (unsigned channel_idx = 0 ; channel_idx<mMetadataNames.size(); channel_idx++)
        {
            // Work out the mean IC50 and Hill in the dataset, and do simulation with this for now.
            // (This will ensure identical behaviour to when there was just one entry)

            double sum = std::accumulate(IC50s[channel_idx].begin(), IC50s[channel_idx].end(), 0.0);
            double mean_ic50 = sum / IC50s[channel_idx].size();
            sum = std::accumulate(hills[channel_idx].begin(), hills[channel_idx].end(), 0.0);
            double mean_hill = sum / hills[channel_idx].size();

            ApplyDrugBlock(mpModel, channel_idx, default_conductances[channel_idx],
                           mConcs[conc_index], mean_ic50, mean_hill);
        }

        double apd90, apd50, upstroke, peak;
        OdeSolution solution = SteadyStatePacingExperiment(mpModel, apd90, apd50, upstroke, peak, 0.1 /*ms printing timestep*/, mConcs[conc_index]);

        // Store some things as member variables for returning later (mainly for testing)
        mApd90s.push_back(apd90); // This is used by TorsadePredict and following method for control.

        InterpolateFromLookupTableForThisConcentration(conc_index);

        if (!DidErrorOccur())
        {
            // Record the control APD90 if this concentration is zero.
            if (fabs(mConcs[conc_index])<1e-12)
            {
                control_apd90 = apd90;
            }
            double delta_apd90 = 100*(apd90 - control_apd90)/control_apd90;
            double lower_delta_apd90, upper_delta_apd90;
            if (mLookupTableAvailable)
            {
                lower_delta_apd90 = 100*(mApd90CredibleRegions[conc_index].first - control_apd90)/control_apd90;
                upper_delta_apd90 = 100*(mApd90CredibleRegions[conc_index].second - control_apd90)/control_apd90;
            }

            if (!mSuppressOutput)
            {
                std::cout << mHertz << "Hz Upstroke velocity = " << upstroke << ", Peak mV = " << peak << ", APD50 = " << apd50 << ", APD90 = " << apd90 << ", percent change APD90 = ";
                if (mLookupTableAvailable)
                {
                    std::cout << lower_delta_apd90 << "," << delta_apd90 << "," << upper_delta_apd90 << "\n";// << std::flush;
                }
                else
                {
                    std::cout << delta_apd90 << "\n";// << std::flush;
                }
            }
            *steady_voltage_results_file_html << "<tr><td>"<< mConcs[conc_index] << "</td><td>" << upstroke << "</td><td>" << peak << "</td><td>" << apd50 << "</td><td>" << apd90 << "</td><td>" << delta_apd90 << "</td></tr>\n";
            *steady_voltage_results_file << mConcs[conc_index] << "\t" << upstroke << "\t" << peak << "\t" << apd50  << "\t" << apd90 << "\t";
            if (mLookupTableAvailable)
            {
                *steady_voltage_results_file << lower_delta_apd90 << "," << delta_apd90 << "," << upper_delta_apd90 << "\n";// << std::flush;
            }
            else
            {
                *steady_voltage_results_file << delta_apd90 << "\n";// << std::flush;
            }
        }
        else
        {
            std::string error_code = GetErrorMessage();
            if (!mSuppressOutput) std::cout << mHertz << "Hz Upstroke velocity = " << error_code << ", Peak mV = " << error_code << ", APD50 = " << error_code << ", APD90 = " << error_code << ", percent change APD90 = " << error_code << "\n";// << std::flush;
            *steady_voltage_results_file_html << "<tr><td>"<< mConcs[conc_index] << "</td><td>" << error_code << "</td><td>" << error_code << "</td><td>" << error_code << "</td><td>" << error_code << "</td><td>" << error_code << "</td></tr>\n";
            *steady_voltage_results_file << mConcs[conc_index] << "\t" << error_code << "\t" << error_code << "\t" << error_code  << "\t" << error_code << "\t";
            if (mLookupTableAvailable)
            {
                *steady_voltage_results_file << error_code << "," << error_code << "," << error_code << "\n";
            }
            else
            {
                *steady_voltage_results_file << error_code << "\n";
            }
        }

        // Create unique filename and write the voltage trace to file...
        std::stringstream filename;
        filename << "conc_" << mConcs[conc_index] << "_voltage_trace.dat";
        boost::shared_ptr<RegularStimulus> p_default_stimulus =
                             boost::static_pointer_cast<RegularStimulus>(mpModel->GetStimulusFunction());
        double s1_period = p_default_stimulus->GetPeriod();
        double s_start = p_default_stimulus->GetStartTime();
        std::vector<double> voltages = solution.GetVariableAtIndex(mpModel->GetSystemInformation()->GetStateVariableIndex("membrane_voltage")); // Voltage should always be zero
        ActionPotentialDownsampler(foldername, filename.str(), solution.rGetTimes(), voltages, s1_period, s_start);
    }// Conc

    // Tidy up
    progress_reporter.PrintFinalising();
    *steady_voltage_results_file_html << "</table>\n</body>\n</html>\n";
    steady_voltage_results_file_html->close();
    steady_voltage_results_file->close();

    mComplete = true;
}

void ApPredictMethods::WriteMessageToFile(const std::string& rMessage)
{
    AbstractActionPotentialMethod::WriteMessageToFile(rMessage);
    assert(mpFileHandler);
    // work out what the absolute path of the message file will be
    const std::string message_file_name = mpFileHandler->GetOutputDirectoryFullPath() + "messages.txt";

    // Check to see if this is the first log message
    FileFinder file_finder(message_file_name, RelativeTo::Absolute);
    bool first_message = !file_finder.Exists();

    // Open the messages file, in appending mode.
    out_stream messages_file = mpFileHandler->OpenOutputFile("messages.txt", std::ios::out | std::ios::app);
    if (first_message)
    {
        *messages_file << "Action potential prediction simulation recorded the following notes:\n";
    }
    *messages_file << " * " << rMessage << "\n";
    messages_file->close();
}

std::vector<double> ApPredictMethods::GetConcentrations(void)
{
    if (!mComplete)
    {
        EXCEPTION("Simulation has not been run - check arguments.");
    }
    return mConcs;
}

std::vector<double> ApPredictMethods::GetApd90s(void)
{
    if (!mComplete)
    {
        EXCEPTION("Simulation has not been run - check arguments.");
    }
    return mApd90s;
}

std::vector<std::pair<double,double> > ApPredictMethods::GetApd90CredibleRegions(void)
{
    if (!mComplete )
    {
        EXCEPTION("Simulation has not been run - check arguments.");
    }

    if (!mLookupTableAvailable)
    {
        EXCEPTION("There was no Lookup Table available for credible interval calculations with these settings.");
    }

    return mApd90CredibleRegions;
}

void ApPredictMethods::ParameterWrapper(boost::shared_ptr<AbstractCvodeCell> pModel, std::vector<std::string>& rMetadataNames)
{
    for (unsigned channel_idx=0; channel_idx<rMetadataNames.size(); channel_idx++)
    {
        if (!pModel->HasParameter(rMetadataNames[channel_idx]))
        {
            WARNING(pModel->GetSystemName() << " does not have '" << rMetadataNames[channel_idx] << "' labelled, please tag it in the CellML file if it is present.");

            // Not all of the models have a distinct fast I_to component.
            // In this case we look for the complete I_to current instead.
            if (rMetadataNames[channel_idx]=="membrane_fast_transient_outward_current_conductance" &&
                pModel->HasParameter("membrane_transient_outward_current_conductance")     )
            {
                WARNING(pModel->GetSystemName() << " does not have 'membrane_fast_transient_outward_current_conductance' labelled, using combined Ito (fast and slow) instead...");
                rMetadataNames[channel_idx] = "membrane_transient_outward_current_conductance";
            }
        }
    }
}



