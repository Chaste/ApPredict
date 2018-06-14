/*

Copyright (c) 2005-2018, University of Oxford.
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

// ApPredict includes
#include "AbstractDataStructure.hpp"
#include "ActionPotentialDownsampler.hpp"
#include "ApPredictMethods.hpp"
#include "BayesianInferer.hpp"
#include "CipaQNetCalculator.hpp"
#include "DoseCalculator.hpp"
#include "LookupTableLoader.hpp"

// Chaste source includes
#include "CommandLineArguments.hpp"
#include "Exception.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "ProgressReporter.hpp"
#include "RegularStimulus.hpp"
#include "SetupModel.hpp"
#include "SteadyStateRunner.hpp"
#include "Timer.hpp"
#include "Warnings.hpp"
#include "ZeroStimulus.hpp"

/**
 * A little helper method that can float around here for now.
 *
 * @param rVec  The vector to return the median of.
 */
double MedianOfStdVectorDouble(const std::vector<double>& rVec)
{
    assert(!rVec.empty());

    std::vector<double> vec = rVec;

    std::sort(vec.begin(), vec.end());

    if (vec.size() % 2u == 0)
    {
        return (vec[vec.size() / 2u - 1u] + vec[vec.size() / 2u]) / 2.0;
    }
    else
    {
        return vec[vec.size() / 2u];
    }
}

/* Citations */
#include "Citations.hpp"

static PetscBool TorsadeCite = PETSC_FALSE;
const char TorsadeCitation[] = "@article{mirams2011simulation,\n"
                               "  title={Simulation of multiple ion channel block provides improved early "
                               "prediction of compounds' clinical torsadogenic risk},\n"
                               "  author={Mirams, G.R. and Cui, Y. and Sher, A. and Fink, M. and "
                               "Cooper, J. and Heath, B.M. and McMahon, N.C. and Gavaghan, D.J. and "
                               "Noble, D.},\n"
                               "  journal={Cardiovascular Research},\n"
                               "  volume={91},\n"
                               "  number={1},\n"
                               "  pages={53--61},\n"
                               "  year={2011},\n"
                               "  doi={10.1093/cvr/CVR044},\n"
                               "}\n";

static PetscBool ApPredictCite = PETSC_FALSE;
const char ApPredictCitation[] = "@article{Williams2015,\n"
                                 "  author = {Williams, Geoff and Mirams, Gary R},\n"
                                 "  doi = {10.1016/j.vascn.2015.05.002},\n"
                                 "  journal = {Journal of pharmacological and toxicological methods},\n"
                                 "  pages = {10--6},\n"
                                 "  title = {A web portal for in-silico action potential predictions},\n"
                                 "  volume = {75},\n"
                                 "  year = {2015},\n"
                                 "}\n";

std::string ApPredictMethods::PrintArguments()
{
    std::string message = "\n**********************************************************************"
                          "*************************\n"
                          "* ApPredict::Please provide some of these inputs:\n*\n"
        + SetupModel::PrintArguments();
    message += PrintCommonArguments();
    return message;
}

std::string ApPredictMethods::PrintCommonArguments()
{
    std::string message = "*\n"
                          "* SPECIFYING PACING:\n"
                          "* --pacing-freq            Pacing frequency (Hz) (optional - defaults "
                          "to 1Hz)\n"
                          "* --pacing-max-time        Maximum time for which to pace the cell "
                          "model in MINUTES\n"
                          "*                          (optional - defaults to time for 10,000 "
                          "paces at this frequency)\n" // Set in AbstractSteadyStateRunner
                          // constructor!
                          "* --pacing-stim-duration   Duration of the square wave stimulus pulse "
                          "applied (ms)\n"
                          "*                          (optional - defaults to stimulus duration "
                          "from CellML)\n"
                          "* --pacing-stim-magnitude  Height of the square wave stimulus pulse "
                          "applied (uA/cm^2)\n"
                          "*                          (optional - defaults to stimulus magnitude "
                          "from CellML)\n"
                          "*\n"
                          "* SPECIFYING DRUG PROPERTIES dose-response properties for each "
                          "channel:\n"
                          "* Channels are named:\n"
                          "* * herg (IKr current - hERG),\n"
                          "* * na (fast sodium current - NaV1.5),\n"
                          "* * nal (late/persistent sodium current - NaV1.5 (perhaps!)),\n"
                          "* * cal (L-type calcium current- CaV1.2),\n"
                          "* * iks (IKs current - KCNQ1 + MinK),\n"
                          "* * ik1 (IK1 current - KCNN4 a.k.a. KCa3.1),\n"
                          "* * ito ([fast] Ito current - Kv4.3 + KChIP2.2).\n"
                          "*\n"
                          "* For each channel you specify dose-response parameters [multiple "
                          "entries for repeat experiments]\n"
                          "*   EITHER with IC50 values (in uM), for example for 'hERG':\n"
                          "* --ic50-herg     hERG IC50    (optional - defaults to \"no effect\")\n"
                          "*   OR with pIC50 values (in log M):\n"
                          "* --pic50-herg    hERG pIC50   (optional - defaults to \"no effect\")\n"
                          "*     (you can use a mixture of these for different channels if you "
                          "wish, \n"
                          "*     e.g. --ic50-herg 16600 --pic50-na 5.3 )\n"
                          "*   AND specify Hill coefficients (dimensionless):\n"
                          "* --hill-herg     hERG Hill    (optional - defaults to \"1.0\")\n"
                          "*   AND specify the saturation effect of the drug on peak conductance "
                          "(%):\n"
                          "* --saturation-herg   saturation level effect of drug (optional - "
                          "defaults to 0%)\n"
                          "*\n"
                          "* SPECIFYING CONCENTRATIONS AT COMMAND LINE:\n"
                          "* --plasma-concs  A list of (space separated) plasma concentrations at "
                          "which to test (uM)\n"
                          "* OR alternatively:\n"
                          "* --plasma-conc-high  Highest plasma concentration to test (uM)\n"
                          "* --plasma-conc-low   Lowest  plasma concentration to test (uM) \n"
                          "*                     (optional - defaults to 0)\n"
                          "*\n"
                          "* both ways of specifying test concentrations have the following "
                          "optional arguments\n"
                          "* --plasma-conc-count  Number of intermediate plasma concentrations to "
                          "test \n"
                          "*                 (optional - defaults to 0 (for --plasma-concs) or 11 "
                          "(for --plasma-conc-high))\n"
                          "* --plasma-conc-logscale <True/False> Whether to use log spacing for "
                          "the plasma concentrations \n"
                          "*\n"
                          "* SPECIFYING CONCENTRATIONS IN A FILE (for PKPD runs):\n"
                          "* if you want to run at concentrations in a file instead of specifying "
                          "at command line, you can do:\n"
                          "* --pkpd-file <relative or absolute filepath>\n"
                          "*   To evaluate APD90s throughout a PKPD profile please provide a file "
                          "with the data format:\n"
                          "*   Time(any "
                          "units)<tab>Conc_trace_1(uM)<tab>Conc_trace_2(uM)<tab>...Conc_trace_N(uM)"
                          "\n"
                          "*   on each row.\n"
                          "*\n"
                          "* UNCERTAINTY QUANTIFICATION:\n"
                          "* --credible-intervals [x y z...] This flag must be present to do "
                          "uncertainty "
                          "calculations. It can optionally be followed by a specific list of "
                          "percentiles that are required\n"
                          "*   (not including 0 or 100, defaults to 95).\n"
                          "* Then to specify 'spread' parameters for assay variability - for use "
                          "with Lookup Tables:\n"
                          "* --pic50-spread-herg      (for each channel that you are providing "
                          "ic50/pic50 values for,\n"
                          "* --hill-spread-herg        herg is just given as an example)\n"
                          "*   (for details of what these spread parameters are see 'sigma' and "
                          "'1/beta' in Table 1 of:\n"
                          "*    Elkins et al. 2013  Journal of Pharmacological and Toxicological \n"
                          "*    Methods, 68(1), 112-122. doi: 10.1016/j.vascn.2013.04.007 )\n"
                          "*\n"
                          "* OTHER OPTIONS:\n"
                          "* --no-downsampling  By default, we print downsampled output to create "
                          "small action potential\n"
                          "*                    traces, but you can switch this off by calling "
                          "this option.\n"
                          "*\n";
    return message;
}

void ApPredictMethods::ReadInIC50HillAndSaturation(
    std::vector<double>& rIc50s, std::vector<double>& rHills,
    std::vector<double>& rSaturations, const unsigned channelIdx)
{
    const std::string channel = mShortNames[channelIdx];
    CommandLineArguments* p_args = CommandLineArguments::Instance();
    bool read_ic50s = false;
    bool read_hills = false;
    bool read_saturations = false;

    // Try loading any arguments given as IC50s
    if (p_args->OptionExists("--ic50-" + channel))
    {
        rIc50s = p_args->GetDoublesCorrespondingToOption("--ic50-" + channel);
        if (p_args->OptionExists("--pic50-" + channel))
        {
            EXCEPTION(
                "Duplicate arguments, you cannot specify both IC50 and pIC50 for "
                << channel << " channel.");
        }
        read_ic50s = true;
    }
    // If those don't exist try loading pIC50s.
    else if (p_args->OptionExists("--pic50-" + channel))
    {
        rIc50s.clear();
        std::vector<double> pIC50s = p_args->GetDoublesCorrespondingToOption("--pic50-" + channel);
        for (unsigned i = 0; i < pIC50s.size(); i++)
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
        if (!(rHills.size() == rIc50s.size()))
        {
            EXCEPTION(
                "If you enter Hill coefficients, there must be one corresponding to "
                "each [p]IC50 measurement.");
        }
        read_hills = true;
    }

    // Try loading any saturations
    if (p_args->OptionExists("--saturation-" + channel))
    {
        rSaturations = p_args->GetDoublesCorrespondingToOption("--saturation-" + channel);
        // But these must correspond to IC50s.
        if (!(rSaturations.size() == rIc50s.size()))
        {
            EXCEPTION(
                "If you enter Saturation levels, there must be one corresponding to "
                "each [p]IC50 measurement.");
        }

        if (rSaturations.size() > 1u)
        {
            WARNING(
                "We haven't yet coded up inference with multiple saturation levels, "
                "just going to use the median value.");
        }
        read_saturations = true;
    }

    // Collect any spread parameter information that has been inputted.
    if (p_args->OptionExists("--pic50-spread-" + channel))
    {
        mPic50Spreads[channelIdx] = p_args->GetDoubleCorrespondingToOption("--pic50-spread-" + channel);
    }
    if (p_args->OptionExists("--hill-spread-" + channel))
    {
        mHillSpreads[channelIdx] = p_args->GetDoubleCorrespondingToOption("--hill-spread-" + channel);
    }
    if (p_args->OptionExists("--saturation-spread-" + channel))
    {
        EXCEPTION(
            "Haven't yet coded up a method to deal with the spread of values on "
            "saturation levels.");
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
        for (unsigned i = 0; i < rIc50s.size(); i++)
        {
            std::cout << rIc50s[i] << " ";
        }
        std::cout << " uM, ";
        if (read_hills)
        {
            std::cout << "Hills = ";
            for (unsigned i = 0; i < rHills.size(); i++)
            {
                std::cout << rHills[i] << ", ";
            }
        }
        else
        {
            std::cout << "Hills = 1.0 (default), ";
        }
        if (read_saturations)
        {
            std::cout << "Saturation levels = ";
            for (unsigned i = 0; i < rSaturations.size(); i++)
            {
                std::cout << rSaturations[i] << " ";
            }
            std::cout << " %." << std::endl;
        }
        else
        {
            std::cout << "Saturation level = 0% (default)." << std::endl;
        }
    }
    else
    {
        std::cout << ": no drug effect\n";
    }
}

void ApPredictMethods::ApplyDrugBlock(
    boost::shared_ptr<AbstractCvodeCell> pModel, unsigned channel_index,
    const double default_conductance, const double concentration,
    const double iC50, const double hill, const double saturation)
{
    // Here we calculate the proportion of the different channels which are still
    // active
    // (at this concentration of this drug)
    const double conductance_factor = AbstractDataStructure::CalculateConductanceFactor(concentration, iC50,
                                                                                        hill, saturation);

    // Some screen output for info.
    if (!mSuppressOutput)
        std::cout << "g_" << mShortNames[channel_index]
                  << " factor = " << conductance_factor
                  << std::endl; // << std::flush;

    // Check the model has this parameter before we try and set it.
    if (pModel->HasParameter(mMetadataNames[channel_index]))
    {
        pModel->SetParameter(mMetadataNames[channel_index],
                             default_conductance * conductance_factor);
    }
    else // We haven't got that conductance parameter, or at least it isn't
    // labelled.
    {
        // If we aren't trying to change it - don't worry, just carry on.
        if (conductance_factor < 1)
        {
            // If the model hasn't got this channel conductance labelled,
            // (but we are trying to change it) throw an error.
            EXCEPTION(
                pModel->GetSystemName()
                << " does not have the current \"" << mMetadataNames[channel_index]
                << "\" labelled, but you have requested a block on this channel.");
        }
    }
}

ApPredictMethods::ApPredictMethods()
        : AbstractActionPotentialMethod(),
          mLookupTableAvailable(false),
          mPercentiles(std::vector<double>{ 2.5, 97.5 }),
          mConcentrationsFromFile(false),
          mComplete(false)
{
    // Here we list the possible drug blocks that can be applied with ApPredict
    mMetadataNames.push_back("membrane_fast_sodium_current_conductance");
    mShortNames.push_back("na");

    mMetadataNames.push_back("membrane_L_type_calcium_current_conductance");
    mShortNames.push_back("cal");

    mMetadataNames.push_back(
        "membrane_rapid_delayed_rectifier_potassium_current_conductance");
    mShortNames.push_back("herg");

    mMetadataNames.push_back(
        "membrane_slow_delayed_rectifier_potassium_current_conductance");
    mShortNames.push_back("iks");

    mMetadataNames.push_back(
        "membrane_inward_rectifier_potassium_current_conductance");
    mShortNames.push_back("ik1");

    mMetadataNames.push_back(
        "membrane_fast_transient_outward_current_conductance");
    mShortNames.push_back("ito");

    mMetadataNames.push_back("membrane_persistent_sodium_current_conductance");
    mShortNames.push_back("nal");

    // We'll use DOUBLE_UNSET to start with for these spread parameters.
    for (unsigned i = 0; i < mMetadataNames.size(); i++)
    {
        mPic50Spreads.push_back(DOUBLE_UNSET);
        mHillSpreads.push_back(DOUBLE_UNSET);
    }

    // There must be a 1:1 mapping between these...
    assert(mMetadataNames.size() == mShortNames.size());

    // Add the fact we're using this code to the citations register
    Citations::Register(TorsadeCitation, &TorsadeCite);
    Citations::Register(ApPredictCitation, &ApPredictCite);

    mProgramName = "Action Potential PreDiCT";
    mOutputFolder = "ApPredict_output/";
}

void ApPredictMethods::SetUpLookupTables()
{
    CommandLineArguments* p_args = CommandLineArguments::Instance();

    if (!p_args->OptionExists("--credible-intervals"))
    {
        // The flag mLookupTableAvailable remains false, and we carry on as normal.
        return;
    }
    else
    {
        if (p_args->GetNumberOfArgumentsForOption("--credible-intervals") > 0)
        {
            // Get list of percentiles to use.
            std::vector<double> percentile_ranges = p_args->GetDoublesCorrespondingToOption("--credible-intervals");
            mPercentiles.clear();
            for (unsigned i = 0; i < percentile_ranges.size(); i++)
            {
                if (percentile_ranges[i] <= 0 || percentile_ranges[i] >= 100)
                {
                    EXCEPTION(
                        "'--credible-intervals' arguments should be given as widths of "
                        "credible interval in percentages. For instance an argument of "
                        "'--credible-intervals 90' will result in 5th and 95th "
                        "percentiles being reported. You specified '"
                        << percentile_ranges[i] << "%' but this number should be more "
                                                   "than zero and less than 100.");
                }
                double remainder_in_tails = 100 - percentile_ranges[i];
                mPercentiles.push_back(0.5 * remainder_in_tails);
                mPercentiles.push_back(100 - 0.5 * remainder_in_tails);
            }
            std::sort(mPercentiles.begin(), mPercentiles.end());
        }
    }

    LookupTableLoader lookup_loader(mpModel->GetSystemName(), this->mHertz);
    if (lookup_loader.IsLookupTableAvailable())
    {
        mpLookupTable = lookup_loader.GetLookupTable();
        mLookupTableAvailable = true;
    }
    else
    {
        WARNING("You asked for '--credible-intervals' but no lookup table is available. Continuing without...");
        mLookupTableAvailable = false;
    }
}

void ApPredictMethods::CalculateDoseResponseParameterSamples(
    const std::vector<std::vector<double> >& rIC50s,
    const std::vector<std::vector<double> >& rHills)
{
    if (!mLookupTableAvailable)
    {
        return;
    }

    /*
       *Prepare an inferred set of IC50s and Hill coefficients
       *for use with the Lookup Table and credible interval calculations.
       */
    mSampledIc50s.resize(mMetadataNames.size());
    mSampledHills.resize(mMetadataNames.size());

    const unsigned num_samples = 1000u;

    // Work out vectors of inferred IC50 and Hills
    // Apply drug block on each channel
    for (unsigned channel_idx = 0; channel_idx < mMetadataNames.size();
         channel_idx++)
    {
        // First just decide whether there is 'no effect here'.
        assert(rIC50s[channel_idx].size() >= 1u);
        if (rIC50s[channel_idx].size() == 1 && rIC50s[channel_idx][0] == -1)
        {
            // No effect here, so just map that out to all random runs
            for (unsigned i = 0; i < num_samples; i++)
            {
                mSampledIc50s[channel_idx].push_back(-1.0);
                mSampledHills[channel_idx].push_back(-1.0);
            }
            // To next channel
            continue;
        }

        // Convert data to pIC50s
        std::vector<double> pIC50s;
        for (unsigned i = 0; i < rIC50s[channel_idx].size(); i++)
        {
            pIC50s.push_back(
                AbstractDataStructure::ConvertIc50ToPic50(rIC50s[channel_idx][i]));
        }

        // Retrieve the Pic50 spread from the command line arguments.
        if (mPic50Spreads[channel_idx] == DOUBLE_UNSET)
        {
            EXCEPTION("No argument --pic50-spread-"
                      << mShortNames[channel_idx]
                      << " has been provided. Cannot calculate credible intervals "
                         "without this.");
        }

        // Infer pIC50 spread.
        BayesianInferer ic50_inferer(PIC50);
        ic50_inferer.SetObservedData(pIC50s);
        ic50_inferer.SetSpreadOfUnderlyingDistribution(mPic50Spreads[channel_idx]);
        ic50_inferer.PerformInference();

        std::vector<double> inferred_pic50s = ic50_inferer.GetSampleMedianValue(
            num_samples); // Get 1000 inferred pIC50s
        for (unsigned i = 0; i < num_samples; i++)
        {
            // Convert pIC50 back to IC50 and store it.
            mSampledIc50s[channel_idx].push_back(
                AbstractDataStructure::ConvertPic50ToIc50(inferred_pic50s[i]));
        }

        // If all Hill coefficient entries are positive, then we will use those for
        // samples too.
        // This long-winded bit of code is just counting how many are positive (must
        // be a better way!).
        bool all_hills_positive = false;
        unsigned temp_counter = 0u;
        for (unsigned i = 0; i < rHills[channel_idx].size(); i++)
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
                WARN_ONCE_ONLY("No argument --hill-spread-"
                               << mShortNames[channel_idx]
                               << " has been provided. "
                                  "Approximating credible intervals without Hill "
                                  "spread info, but you will get better answers with "
                                  "it.");

                // If we can't guess the hill spread, then just use mean Hill that we
                // have.
                std::vector<double> hills_this_channel = rHills[channel_idx];
                double mean_hill = std::accumulate(hills_this_channel.begin(),
                                                   hills_this_channel.end(), 0.0)
                    / hills_this_channel.size();

                for (unsigned i = 0; i < num_samples; i++)
                {
                    // We aren't going to attempt to do inference on Hills, just IC50s,
                    // push back mean Hill.
                    mSampledHills[channel_idx].push_back(mean_hill);
                }
            }
            else
            {
                // Infer Hill spread.
                BayesianInferer hill_inferer(HILL);
                hill_inferer.SetObservedData(rHills[channel_idx]);
                // This works with the Beta parameter, not the 1/Beta. So do 1/1/Beta to
                // get Beta back!
                hill_inferer.SetSpreadOfUnderlyingDistribution(
                    1.0 / mHillSpreads[channel_idx]);
                hill_inferer.PerformInference();

                mSampledHills[channel_idx] = hill_inferer.GetSampleMedianValue(
                    num_samples); // Get 1000 inferred Hills
            }
        }
        else
        {
            // There have been not enough 'real' Hill coefficients specified,
            // so assume they are all missing (otherwise inference would go mad
            // with some positive and some negative).
            for (unsigned i = 0; i < num_samples; i++)
            {
                // We aren't going to attempt to do inference on Hills,
                // just push back 'no measurement' for now.
                mSampledHills[channel_idx].push_back(-1.0);
            }
        }
    }
}

void ApPredictMethods::InterpolateFromLookupTableForThisConcentration(
    const unsigned concIndex,
    const std::vector<double>& rMedianSaturationLevels)
{
    // If we don't have a lookup table, we aren't going to do confidence
    // intervals.
    if (!mLookupTableAvailable)
    {
        return;
    }

    const unsigned table_dim = mpLookupTable->GetDimension();

    std::vector<double> credible_intervals;

    // If this is the first concentration (control) say the percent change must be
    // zero
    // or otherwise a small interpolation error will result
    // (from potentially running to different steady state with --pacing-max-time
    // <x> ).
    if (concIndex == 0u)
    {
        for (unsigned i = 0; i < mPercentiles.size(); i++)
        {
            credible_intervals.push_back(mApd90s[concIndex]);
        }
        mApd90CredibleRegions[concIndex] = credible_intervals;
        return;
    }

    // The first channel entry in mSampledIc50s will give us the number of random
    // samples.
    const unsigned num_samples = mSampledIc50s[0].size();

    // In the lookup table the order of parameters is given in the filename:
    // "4d_hERG_IKs_INa_ICaL_generator.arch"
    std::vector<std::string> required_channels;
    required_channels.push_back("herg");
    if (table_dim > 3u)
    {
        required_channels.push_back("iks");
    }
    required_channels.push_back("na");
    required_channels.push_back("cal");

    // This slightly complicated loop is just seeing which entry in
    // mSampledIC50/Hills corresponds
    // to the ones that we want, so we've listed the ones we want above and search
    // for them in
    // mShortNames here.
    std::map<unsigned, unsigned> map_to_metadata_idx;
    for (unsigned channel_idx = 0; channel_idx < table_dim; channel_idx++)
    {
        for (unsigned i = 0; i < mShortNames.size(); i++)
        {
            if (mShortNames[i] == required_channels[channel_idx])
            {
                map_to_metadata_idx[channel_idx] = i;
                break;
            }
        }
    }
    assert(map_to_metadata_idx.size() == table_dim);

    std::cout << "Calculating confidence intervals from Lookup Table...";
    std::vector<std::vector<double> > sampling_points;
    for (unsigned rand_idx = 0; rand_idx < num_samples; rand_idx++)
    {
        std::vector<double> sample_required_at(table_dim);
        for (unsigned i = 0; i < table_dim; i++)
        {
            sample_required_at[i] = AbstractDataStructure::CalculateConductanceFactor(
                mConcs[concIndex], mSampledIc50s[map_to_metadata_idx[i]][rand_idx],
                mSampledHills[map_to_metadata_idx[i]][rand_idx],
                rMedianSaturationLevels[map_to_metadata_idx[i]]);
        }
        sampling_points.push_back(sample_required_at);
    }

    std::vector<std::vector<double> > predictions = mpLookupTable->Interpolate(sampling_points);
    assert(predictions.size() == mSampledIc50s[0].size());

    /*
       * Compile all the lookup table predictions into a vector that we can sort
   * to
       * get percentiles.
       */
    // std::vector<double> apd_50_predictions;
    std::vector<double> apd_90_predictions;
    for (unsigned rand_idx = 0; rand_idx < num_samples; rand_idx++)
    {
        // apd_50_predictions.push_back(predictions[rand_idx][1]);
        apd_90_predictions.push_back(predictions[rand_idx][0]);
    }
    std::sort(apd_90_predictions.begin(), apd_90_predictions.end());

    for (unsigned i = 0; i < mPercentiles.size(); i++)
    {
        // Now work out the confidence intervals, err on conservative side.
        unsigned index_in_sorted_apd90_vector;
        if (mPercentiles[i] < 50)
        {
            index_in_sorted_apd90_vector = floor(mPercentiles[i] / 100.0 * (double)(num_samples));
        }
        else
        {
            index_in_sorted_apd90_vector = ceil(mPercentiles[i] / 100.0 * (double)(num_samples));
        }
        credible_intervals.push_back(
            apd_90_predictions[index_in_sorted_apd90_vector]);
    }

    mApd90CredibleRegions[concIndex] = credible_intervals;
    std::cout << "done." << std::endl;
}

void ApPredictMethods::Run()
{
    // Make and clean the above directories.
    mpFileHandler.reset(new OutputFileHandler(mOutputFolder));

    // This class will get model definition from command line, so we don't pass in
    // model index.
    SetupModel setup(this->mHertz, UNSIGNED_UNSET, mpFileHandler);
    mpModel = setup.GetModel();

    SetUpLookupTables();

    CommonRunMethod();
}

void ApPredictMethods::CommonRunMethod()
{
    // Arguments that take default values
    std::vector<std::vector<double> > IC50s;
    std::vector<std::vector<double> > hills;
    std::vector<std::vector<double> > saturations;

    std::vector<double> unset;
    unset.push_back(-1); // -1 is our code for 'unset'.

    // Set the initial values of these
    for (unsigned i = 0; i < mMetadataNames.size(); i++)
    {
        IC50s.push_back(unset);
        hills.push_back(unset);
        saturations.push_back(unset);
    }

    if (CommandLineArguments::Instance()->OptionExists("--pkpd-file"))
    {
        FileFinder pkpd_file(
            CommandLineArguments::Instance()->GetStringCorrespondingToOption(
                "--pkpd-file"),
            RelativeTo::AbsoluteOrCwd);
        if (!pkpd_file.IsFile())
        {
            EXCEPTION(
                "The File '"
                << pkpd_file.GetAbsolutePath()
                << "' does not exist. Please give a relative or absolute path.");
        }

        if (CommandLineArguments::Instance()->OptionExists("--plasma-conc-high"))
        {
            EXCEPTION(
                "The argument --plasma-conc-high will be ignored. Using PKPD file to "
                "set concentrations. Please remove it to avoid confusion!");
        }

        if (CommandLineArguments::Instance()->OptionExists("--plasma-concs"))
        {
            EXCEPTION(
                "The argument --plasma-concs will be ignored. Using PKPD file to set "
                "concentrations. Please remove it to avoid confusion!");
        }

        // Set up a structure to read the PK concentrations in.
        mpPkpdReader = boost::shared_ptr<PkpdDataStructure>(new PkpdDataStructure(pkpd_file));
        mConcentrationsFromFile = true;

        // Calculate maximum concentration to use... add 10% to the maximum we saw.
        DoseCalculator dose_calculator(1.1 * mpPkpdReader->GetMaximumConcentration());
        dose_calculator.SetNumSubdivisions(97); // Loads of detail for these sims
        // to be accurately interpolated
        // later.
        mConcs = dose_calculator.GetConcentrations();
    }
    else
    // this default DoseCalculator constructor reads command line arguments
    // to set the plasma concentrations.
    {
        DoseCalculator dose_calculator;
        mConcs = dose_calculator.GetConcentrations();
    }

    // We check the desired parameters are present in the model, warn if not.
    // This method also changes some metadata names if the model has variants that
    // will do, but aren't ideal
    // and warns if it does this.
    ParameterWrapper(mpModel, mMetadataNames);

    // Use a helper method to read in IC50 from either --ic50 or --pic50
    // arguments.
    // Note this is now in micro Molar (1x10^-6 Molar) as per most Pharma use.
    for (unsigned channel_idx = 0; channel_idx < mMetadataNames.size();
         channel_idx++)
    {
        ReadInIC50HillAndSaturation(IC50s[channel_idx], hills[channel_idx],
                                    saturations[channel_idx], channel_idx);
    }

    if (!mSuppressOutput)
    {
        std::cout << "* max free plasma concentration = " << mConcs.back()
                  << " uM\n"
                     "* min free plasma concentration = "
                  << mConcs[0] << " uM\n"
                                  "* number of plasma concentrations = "
                  << mConcs.size() << "\n"; // << std::flush;
    }

    // The following names are fixed and correspond to metadata tags.
    // We record the default parameter values that the model uses.
    // All the drug block models should include these parameter labels
    // But you only get a warning if not, so check the warnings...
    std::vector<double> default_conductances;
    for (unsigned channel_idx = 0; channel_idx < mMetadataNames.size();
         channel_idx++)
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
    boost::shared_ptr<RegularStimulus> p_reg_stim = boost::static_pointer_cast<RegularStimulus>(
        mpModel->GetStimulusFunction());
    p_reg_stim->SetStartTime(5.0);

    // If we are using ORdCiPAv1 and 0.5Hz, calculate qNet.
    bool calculate_qNet = false;
    out_stream q_net_results_file;
    if (model_name == "ohara_rudy_cipa_v1_2017" && p_reg_stim->GetPeriod() == 2000)
    {
        calculate_qNet = true;
        q_net_results_file = mpFileHandler->OpenOutputFile("q_net.txt");
        *q_net_results_file << "Concentration(uM)\tqNet(C/F)" << std::endl;
    }

    // Print out a progress file for monitoring purposes.
    ProgressReporter progress_reporter(mOutputFolder, 0.0,
                                       (double)(mConcs.size()));
    progress_reporter.PrintInitialising();

    // Open files and write headers
    out_stream steady_voltage_results_file_html = mpFileHandler->OpenOutputFile("voltage_results.html");

    out_stream steady_voltage_results_file = mpFileHandler->OpenOutputFile("voltage_results.dat");
    *steady_voltage_results_file << "Concentration(uM)\tUpstrokeVelocity(mV/"
                                    "ms)\tPeakVm(mV)\tAPD50(ms)\tAPD90(ms)\t";
    if (mLookupTableAvailable)
    {
        for (unsigned i = 0; i < mPercentiles.size(); i++)
        {
            std::string lower_or_upper = "low";
            if (mPercentiles[i] > 50)
            {
                lower_or_upper = "upp";
                if (mPercentiles[i - 1] < 50)
                {
                    *steady_voltage_results_file << "median_delta_APD90"
                                                 << ",";
                }
            }
            double credible_interval;
            if (mPercentiles[i] < 50)
            {
                credible_interval = 100 - 2 * mPercentiles[i];
            }
            else
            {
                credible_interval = 100 - 2 * (100 - mPercentiles[i]);
            }
            *steady_voltage_results_file << "dAp" << credible_interval << "%"
                                         << lower_or_upper;
            if (i < mPercentiles.size() - 1u)
            {
                *steady_voltage_results_file << ",";
            }
        }
        *steady_voltage_results_file << std::endl;
    }
    else
    {
        *steady_voltage_results_file << "delta_APD90(%)\n";
    }

    *steady_voltage_results_file_html << "<html>\n<head><title>" << mProgramName
                                      << " results</title></head>\n";
    *steady_voltage_results_file_html
        << "<STYLE TYPE=\"text/css\">\n<!--\nTD{font-size: "
           "12px;}\n--->\n</STYLE>\n";
    *steady_voltage_results_file_html << "<body>\n";
    *steady_voltage_results_file_html
        << "<table width=\"60%\" style=\"background-color:white\" border=\"1\" "
           "cellpadding=\"2\" cellspacing=\"0\">\n";
    *steady_voltage_results_file_html << "<tr><td>Concentration "
                                         "(uM)</td><td>Upstroke Velocity "
                                         "(mV/ms)</td><td>Peak membrane voltage "
                                         "(mV)</td><td>APD50 (ms)</td><td>APD90 "
                                         "(ms)</td><td>Change in APD90 "
                                         "(%)</td></tr>\n"; // Header line

    /*
   * Work out the median IC50, Hill and saturation to use if more than one were
   * provided
   */
    std::vector<double> median_ic50; // vector is over channel indices
    std::vector<double> median_hill; //               ""
    std::vector<double> median_saturation; //               ""
    for (unsigned channel_idx = 0; channel_idx < mMetadataNames.size();
         channel_idx++)
    {
        // If we only have one IC50 and Hill value specified, then just use them and
        // move on.
        if (IC50s[channel_idx].size() == 1u && hills[channel_idx].size() == 1u)
        {
            median_ic50.push_back(IC50s[channel_idx][0]);
            median_hill.push_back(hills[channel_idx][0]);
            median_saturation.push_back(saturations[channel_idx][0]);
        }
        else
        {
            // Otherwise we have two scenarios
            // 1. We know something about the spread and should take the median of the
            // inferred distribution samples
            //    This should ensure that the predicted line is at the 50th percentile
            //    of the credible interval predictions.
            // 2. We know nothing about the spread, and should just take the median of
            // the multiple values.
            if (mLookupTableAvailable)
            {
                // Work out the median inferred IC50 and hill from the multiple values
                // dataset,
                // and do the simulation with those. Note that the median will give the
                // same value
                // whether we use IC50s or pIC50s, whereas the mean is skewed...
                median_ic50.push_back(
                    MedianOfStdVectorDouble(mSampledIc50s[channel_idx]));
                median_hill.push_back(
                    MedianOfStdVectorDouble(mSampledHills[channel_idx]));

                // TODO: Clever way to do inference on saturation levels too. No data
                // analysed to work out
                // their distributions yet, so we'll just use the median for now as
                // below.
            }
            else
            {
                // Work out the median pIC50 and Hill in the multiple values dataset,
                // and do simulation with these.
                // This is our best estimate of the above more accurate version when we
                // have no other information.
                std::vector<double> pIC50s;
                for (unsigned i = 0; i < IC50s[channel_idx].size(); i++)
                {
                    pIC50s.push_back(
                        AbstractDataStructure::ConvertIc50ToPic50(IC50s[channel_idx][i]));
                }
                median_ic50.push_back(AbstractDataStructure::ConvertPic50ToIc50(
                    MedianOfStdVectorDouble(pIC50s)));
                median_hill.push_back(MedianOfStdVectorDouble(hills[channel_idx]));
            }

            // We've no clever way of dealing with this yet, just take median of
            // saturation levels and use that all the time.
            median_saturation.push_back(
                MedianOfStdVectorDouble(saturations[channel_idx]));
        }
    }

    /*
     * START LOOP OVER EACH CONCENTRATION TO TEST WITH
     */
    bool reliable_credible_intervals = true;
    mApd90CredibleRegions.resize(mConcs.size());
    double control_apd90 = 0;
    for (unsigned conc_index = 0u; conc_index < mConcs.size(); conc_index++)
    {
        progress_reporter.Update((double)(conc_index));
        std::cout << "Drug Conc = " << mConcs[conc_index] << " uM"
                  << std::endl; //<< std::flush;

        // Apply drug block on each channel
        for (unsigned channel_idx = 0; channel_idx < mMetadataNames.size();
             channel_idx++)
        {
            ApplyDrugBlock(mpModel, channel_idx, default_conductances[channel_idx],
                           mConcs[conc_index], median_ic50[channel_idx],
                           median_hill[channel_idx], median_saturation[channel_idx]);
        }

        double apd90, apd50, upstroke, peak, peak_time, ca_max, ca_min;
        OdeSolution solution = SteadyStatePacingExperiment(
            mpModel, apd90, apd50, upstroke, peak, peak_time, ca_max, ca_min,
            0.1 /*ms printing timestep*/, mConcs[conc_index]);

        if (calculate_qNet)
        {
            CipaQNetCalculator calculator(mpModel);
            double q_net = calculator.ComputeQNet();
            std::cout << "qNet at " << mConcs[conc_index] << "uM = " << q_net
                      << " C/F" << std::endl;
            *q_net_results_file << mConcs[conc_index] << "\t" << q_net << std::endl;

            if (conc_index == mConcs.size() - 1u && this->GetMaxNumPaces() < 750u)
            {
                std::stringstream message;
                message << "Warning: qNet is calculated after at least 750 paces in "
                           "FDA publications. You are doing "
                        << this->GetMaxNumPaces()
                        << " paces at " << mConcs[conc_index] << "uM, increase maximum pacing time if using these "
                                                                 "simulation results for CiPA purposes.";
                WriteMessageToFile(message.str());
            }

            if (q_net == std::numeric_limits<double>::quiet_NaN())
            {
                std::stringstream message;
                message << "At a concentration of " << mConcs[conc_index] << "uM qNet was not calculated as the AP did not repolarise (this indicates very high risk).";
                WriteMessageToFile(message.str());
            }
        }

        // Store some things as member variables for returning later (mainly for
        // testing)
        mApd90s.push_back(apd90); // This is used by TorsadePredict and following
        // method for control.

        InterpolateFromLookupTableForThisConcentration(conc_index,
                                                       median_saturation);

        if (!DidErrorOccur())
        {
            // Record the control APD90 if this concentration is zero.
            if (fabs(mConcs[conc_index]) < 1e-12)
            {
                control_apd90 = apd90;
            }
            double delta_apd90 = 100 * (apd90 - control_apd90) / control_apd90;
            std::vector<double> delta_percentiles(mPercentiles.size());
            if (mLookupTableAvailable)
            {
                for (unsigned i = 0; i < mPercentiles.size(); i++)
                {
                    delta_percentiles[i] = 100 * (mApd90CredibleRegions[conc_index][i] - control_apd90) / control_apd90;
                }
            }

            if (!mSuppressOutput)
            {
                std::cout << mHertz << "Hz Upstroke velocity = " << upstroke
                          << ", Peak mV = " << peak << ", APD50 = " << apd50
                          << ", APD90 = " << apd90 << ", percent change APD90 = ";
                if (mLookupTableAvailable)
                {
                    std::cout << delta_percentiles[0] << "," << delta_apd90 << ","
                              << delta_percentiles[mPercentiles.size() - 1u]
                              << std::endl; // << std::flush;
                }
                else
                {
                    std::cout << delta_apd90 << std::endl; // << std::flush;
                }
            }
            *steady_voltage_results_file_html
                << "<tr><td>" << mConcs[conc_index] << "</td><td>" << upstroke
                << "</td><td>" << peak << "</td><td>" << apd50 << "</td><td>" << apd90
                << "</td><td>" << delta_apd90 << "</td></tr>\n";
            *steady_voltage_results_file << mConcs[conc_index] << "\t" << upstroke
                                         << "\t" << peak << "\t" << apd50 << "\t"
                                         << apd90 << "\t";
            if (mLookupTableAvailable)
            {
                for (unsigned i = 0; i < mPercentiles.size(); i++)
                {
                    if (mPercentiles[i] > 50 && mPercentiles[i - 1] < 50)
                    {
                        *steady_voltage_results_file << delta_apd90 << ",";
                    }
                    // Now add a check to see whether the middle 30% of our credible interval contains
                    // the simulated 'median' answer. If not, log it and report a warning later.
                    if ((mPercentiles[i] < 35 && delta_percentiles[i] > delta_apd90)
                        || (mPercentiles[i] > 75 && delta_percentiles[i] < delta_apd90))
                    {
                        reliable_credible_intervals = false;
                    }
                    *steady_voltage_results_file << delta_percentiles[i];
                    if (i < mPercentiles.size() - 1u)
                    {
                        *steady_voltage_results_file << ",";
                    }
                }
                *steady_voltage_results_file << std::endl; // << std::flush;
            }
            else
            {
                *steady_voltage_results_file << delta_apd90
                                             << std::endl; // << std::flush;
            }
        }
        else
        {
            std::string error_code = GetErrorMessage();
            if (!mSuppressOutput)
                std::cout << mHertz << "Hz Upstroke velocity = " << error_code
                          << ", Peak mV = " << error_code << ", APD50 = " << error_code
                          << ", APD90 = " << error_code
                          << ", percent change APD90 = " << error_code
                          << "\n"; // << std::flush;
            *steady_voltage_results_file_html
                << "<tr><td>" << mConcs[conc_index] << "</td><td>" << error_code
                << "</td><td>" << error_code << "</td><td>" << error_code
                << "</td><td>" << error_code << "</td><td>" << error_code
                << "</td></tr>\n";
            *steady_voltage_results_file << mConcs[conc_index] << "\t" << error_code
                                         << "\t" << error_code << "\t" << error_code
                                         << "\t" << error_code << "\t";
            if (mLookupTableAvailable)
            {
                *steady_voltage_results_file << error_code << "," << error_code << ","
                                             << error_code << std::endl;
            }
            else
            {
                *steady_voltage_results_file << error_code << std::endl;
            }
        }

        // Create unique filename and write the voltage trace to file...
        std::stringstream filename;
        filename << "conc_" << mConcs[conc_index] << "_voltage_trace.dat";
        boost::shared_ptr<RegularStimulus> p_default_stimulus = boost::static_pointer_cast<RegularStimulus>(mpModel->GetStimulusFunction());
        double s1_period = p_default_stimulus->GetPeriod();
        double s_start = p_default_stimulus->GetStartTime();
        std::vector<double> voltages = solution.GetVariableAtIndex(mpModel->GetSystemInformation()->GetStateVariableIndex("membrane_voltage"));
        double window = s1_period;
        if (this->mPeriodTwoBehaviour)
        {
            window *= 2.0;
        }
        ActionPotentialDownsampler(mOutputFolder, filename.str(), solution.rGetTimes(), voltages, window, s_start);
    } // Conc

    if (!reliable_credible_intervals)
    {
        WriteMessageToFile("Warning: the credible intervals here (from lookup tables) do not align with simulation - treat them with caution, and ideally report simulation details to allow us to refine lookup tables.");
    }

    // Tidy up
    progress_reporter.PrintFinalising();
    *steady_voltage_results_file_html << "</table>\n</body>\n</html>\n";
    steady_voltage_results_file_html->close();
    steady_voltage_results_file->close();
    if (calculate_qNet)
    {
        q_net_results_file->close();
    }

    if (mConcentrationsFromFile)
    {
        // Copy the input file to the results folder, for posterity...
        FileFinder pkpd_file(
            CommandLineArguments::Instance()->GetStringCorrespondingToOption(
                "--pkpd-file"),
            RelativeTo::AbsoluteOrCwd);
        mpFileHandler->CopyFileTo(pkpd_file);

        // Open a results file - in a try catch as it is conceivable someone could
        // have named their PK file this!
        out_stream p_output_file;
        try
        {
            p_output_file = mpFileHandler->OpenOutputFile("pkpd_results.txt");
        }
        catch (Exception& e)
        {
            EXCEPTION("ApPredict could not open a new output file called pkpd_results.txt. Error was: '" << e.GetMessage() << "'");
        }
        *p_output_file << "Time";
        for (unsigned i = 0; i < mpPkpdReader->GetNumberOfPatients(); i++)
        {
            *p_output_file << "\tAPD90_for_patient_" << i << "(ms)";
        }
        *p_output_file << std::endl;

        std::vector<std::string> times = mpPkpdReader->GetTimes();
        for (unsigned i = 0; i < times.size(); i++)
        {
            *p_output_file << times[i];
            const std::vector<double>& r_concs_at_this_time = mpPkpdReader->GetConcentrationsAtTimeIndex(i);
            for (unsigned p = 0; p < r_concs_at_this_time.size(); p++)
            {
                double interpolated_apd90 = DoLinearInterpolation(r_concs_at_this_time[p], mConcs, mApd90s);
                *p_output_file << "\t" << interpolated_apd90;
            }
            *p_output_file << std::endl;
        }
        p_output_file->close();
    }
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
    out_stream messages_file = mpFileHandler->OpenOutputFile(
        "messages.txt", std::ios::out | std::ios::app);
    if (first_message)
    {
        *messages_file << "Action potential prediction simulation recorded the "
                          "following notes:\n";
    }
    *messages_file << " * " << rMessage << std::endl;
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

std::vector<std::vector<double> > ApPredictMethods::GetApd90CredibleRegions(
    void)
{
    if (!mComplete)
    {
        EXCEPTION("Simulation has not been run - check arguments.");
    }

    if (!mLookupTableAvailable)
    {
        EXCEPTION(
            "There was no Lookup Table available for credible interval "
            "calculations with these settings.");
    }

    return mApd90CredibleRegions;
}

void ApPredictMethods::ParameterWrapper(
    boost::shared_ptr<AbstractCvodeCell> pModel,
    std::vector<std::string>& rMetadataNames)
{
    for (unsigned channel_idx = 0; channel_idx < rMetadataNames.size();
         channel_idx++)
    {
        if (!pModel->HasParameter(rMetadataNames[channel_idx]))
        {
            WARNING(
                pModel->GetSystemName()
                << " does not have '" << rMetadataNames[channel_idx]
                << "' labelled, please tag it in the CellML file if it is present.");

            // Not all of the models have a distinct fast I_to component.
            // In this case we look for the complete I_to current instead.
            if (rMetadataNames[channel_idx] == "membrane_fast_transient_outward_current_conductance" && pModel->HasParameter("membrane_transient_outward_current_conductance"))
            {
                WARNING(pModel->GetSystemName()
                        << " does not have "
                           "'membrane_fast_transient_outward_current_conductance' "
                           "labelled, using combined Ito (fast and slow) instead...");
                rMetadataNames[channel_idx] = "membrane_transient_outward_current_conductance";
            }
        }
    }
}

void ApPredictMethods::SetOutputDirectory(const std::string& rOuputDirectory)
{
    mOutputFolder = rOuputDirectory;
}

/* Perform linear interpolation to get an estimate of y_star at x_star */
double ApPredictMethods::DoLinearInterpolation(
    double x_star, const std::vector<double>& rX,
    const std::vector<double>& rY) const
{
    if (x_star <= rX[0])
    {
        return rY[0];
    }
    if (x_star >= rX.back())
    {
        return rY.back();
    }

    auto lower = std::lower_bound(rX.cbegin(), rX.cend(), x_star);
    unsigned lower_idx = lower - rX.cbegin();

    // I'll let the compiler tidy all this up!
    double lower_x = rX[lower_idx - 1u];
    double upper_x = rX[lower_idx];
    double lower_y = rY[lower_idx - 1u];
    double upper_y = rY[lower_idx];

    return lower_y + ((x_star - lower_x) / (upper_x - lower_x)) * (upper_y - lower_y);
}
