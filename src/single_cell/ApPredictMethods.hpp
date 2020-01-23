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

#ifndef APPREDICTMETHODS_HPP_
#define APPREDICTMETHODS_HPP_

#include "AbstractActionPotentialMethod.hpp"
#include "AbstractCvodeCell.hpp"
#include "LookupTableGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "PkpdDataStructure.hpp"

/**
 * Common code to allow this to be run as a test via scons and also as a
 * executable application from the apps folder.
 *
 * This class provides all of the functionality behind the ApPredict executable
 * and is used to generate simple APD90 predictions based on conductance block
 * of the ion channels specified in the constructor.
 *
 * The class can also use a Lookup Table to generate thousands of APD90 predictions
 * very very quickly, to create Bayesian 'credible intervals' (a bit like confidence
 * intervals) around the model APD90 predictions, at each concentration.
 * These are based on the input data uncertainty, as characterised in the paper:
 * Elkins et al. 2013  Journal of Pharmacological and Toxicological
 * Methods, 68(1), 112-122. doi: 10.1016/j.vascn.2013.04.007
 *
 * For most details of using this class, please compile and run the ApPredict executable,
 * and it will provide a 'help' style output of the command line arguments to use.
 */
class ApPredictMethods : public AbstractActionPotentialMethod
{
private:
    friend class TestPkpdInterpolations; // To test linear interpolation private method.
    /**
     * A method to avoid lots of copying and pasting in the main drug block application method.
     *
     * @param pModel  The CVODE model.
     * @param channel_index  The index of the channel (in mMetadataNames) that we are blocking here.
     * @param default_conductance  The default value of the conductance before we started messing.
     * @param concentration  The current drug concentration.
     * @param iC50  The IC50 for this channel.
     * @param hill  The Hill for this channel.
     * @param saturation  The saturation level for this channel and drug (e.g. 0% = full block, 150% = 50% activator).
     */
    void ApplyDrugBlock(boost::shared_ptr<AbstractCvodeCell> pModel,
                        unsigned channel_index,
                        const double default_conductance,
                        const double concentration,
                        const double iC50,
                        const double hill,
                        const double saturation);

    /**
     * Takes the inputted IC50 and Hill coefficients
     * and uses information on the spread of the particular assay
     * to infer a probability distribution for the 'true' underlying
     * median IC50 and Hill.
     *
     * Then stores samples from the inferred PDF in #mSampledIc50s
     * and #mSampledHills.
     *
     * @param rIC50s  The IC50 values for each channel, and any repeated measurements (inner vec).
     * @param rHills  The Hill coefficients for each channel, and any repeated measurements (inner vec).
     */
    void CalculateDoseResponseParameterSamples(const std::vector<std::vector<double> >& rIC50s,
                                               const std::vector<std::vector<double> >& rHills);

    /**
     * Uses the lookup table and entries in mSampledIc50s and mSampledHills
     * to generate a probability distribution of APD90 predictions. This is then
     * stored in mAllApd90s, and mApd90CredibleRegions is populated.
     *
     * @param concIndex  The index of the concentration (in mConcs).
     * @param rMedianSaturationLevels  The saturation levels for each channel to assume in all samples.
     */
    void InterpolateFromLookupTableForThisConcentration(const unsigned concIndex,
                                                        const std::vector<double>& rMedianSaturationLevels);

    /**
      * Perform linear interpolation to get an estimate of y_star at x_star
      * @param x_star The independent variable at which to get an interpolated value
      * @param rX  The vector of independent variables.
      * @param rY  The vector of dependent variables to interpolate between.
      */
    double DoLinearInterpolation(double x_star, const std::vector<double>& rX, const std::vector<double>& rY) const;

    /** The Oxford metadata names of the conductances we may modify with this class */
    std::vector<std::string> mMetadataNames;

    /** Shortened versions of the names, corresponding to the input argument names */
    std::vector<std::string> mShortNames;

    /**
     * The IC50 samples for credible interval calculations.
     * The first index is for channel (corresponding to mMetadataNames)
     * The second (inner) vector is for each random sample.
     */
    std::vector<std::vector<double> > mSampledIc50s;

    /**
     * The Hill coefficient samples for credible interval calculations.
     * The first index is for channel (corresponding to mMetadataNames)
     * The second (inner) vector is for each random sample.
     */
    std::vector<std::vector<double> > mSampledHills;

    /** The inputted spread parameters - for pIC50 Logistic Distbn this is 'sigma' */
    std::vector<double> mPic50Spreads;

    /** The inputted spread parameters - for Hill Log-Logisitic this is '1/Beta' */
    std::vector<double> mHillSpreads;

    /** Whether there is a lookup table we can use for credible interval calculations */
    bool mLookupTableAvailable;

    /** A pointer to a lookup table */
    boost::shared_ptr<AbstractUntemplatedLookupTableGenerator> mpLookupTable;

    /**
     * A vector of pairs used to store the credible regions for APD90s,
     * calculated in the main method if a suitable Lookup Table is present.
     *
     * The outer vector loops over concentrations.
     * At each concentration we have a vector of values for the percentiles in #mPercentiles.
     */
    std::vector<std::vector<double> > mApd90CredibleRegions;
    
    /**
     * The percentiles that the credible region APD90 values in #mApd90CredibleRegions
     * correspond to.
     */
    std::vector<double> mPercentiles;

    /**
     * Whether we are running a Pharmacokinetics simulation with concentrations read from file.
     */
    bool mConcentrationsFromFile;

    /** A data reader class to hold information read from the PKPD file */
    boost::shared_ptr<PkpdDataStructure> mpPkpdReader;

protected:
    /** Whether the simulation completed successfully */
    bool mComplete;

    /** A vector used to store the APD90s calculated in the main method */
    std::vector<double> mApd90s;

    /** A vector used to store the Drug Concentrations at which APDs are calculated*/
    std::vector<double> mConcs;

    /** A file handler that points to the output directory */
    boost::shared_ptr<OutputFileHandler> mpFileHandler;

    /** The program name, for messages, refreshed on each Run call. */
    std::string mProgramName;

    /** The output folder, refreshed on each Run call. */
    std::string mOutputFolder;

    /** The model we're working with, refreshed on each Run call.*/
    boost::shared_ptr<AbstractCvodeCell> mpModel;

    /**
     * Read any input arguments corresponding to a particular channel and calculate the IC50 value (in uM)
     * from either raw IC50 (in uM) or pIC50 (in M).
     *
     * Also stores any 'spread' information (for credible interval calculations) in the mSampledIc50s and
     * mSampledHills member variables.
     *
     * @param rIc50  default IC50 values (usually -1), overwritten if argument is present to a value
     *               (assumed to be in uM).
     * @param rHill  default Hill coefficients (usually -1), overwritten if argument is present.
     * @param rSaturation  default saturation levels (usually -1), overwritten if argument is present.
     * @param channelIdx  The index of the channel in mMetadataNames.
     */
    void ReadInIC50HillAndSaturation(std::vector<double>& rIc50,
                                     std::vector<double>& rHill,
                                     std::vector<double>& rSaturation,
                                     const unsigned channelIdx);

    /**
     * Write a log message to the messages.txt file that should be displayed alongside the results
     * (for example a warning that the cell failed to de/re-polarise at a certain concentration.)
     * @param rMessage  The message to write to file.
     */
    void WriteMessageToFile(const std::string& rMessage);

    /**
     * Converts a set of intended metadata names into the ones we should use...
     *
     * e.g. if the model doesn't have a particular variant of a current, but does have a similar one.
     * Provides warnings if it ever changes the metadata names underneath you.
     *
     * @param pModel  the model
     * @param rMetadataNames  The metadata names we want to use, these may be changed.
     */
    void ParameterWrapper(boost::shared_ptr<AbstractCvodeCell> pModel, std::vector<std::string>& rMetadataNames);

    /**
     * This is a helper method to print out the available arguments.
     *
     * @return the arguments ApPredict, TorsadePredict, PkpdInterpolator all take.
     */
    static std::string PrintCommonArguments();

    /**
     * A run method common to both ApPredict, TorsadePredict and PkpdInterpolator.
     *
     * Does the majority of the work!
     */
    void CommonRunMethod();

    /**
     * A method to look for / download and unarchive any lookup tables.
     *
     * This method contains the following logic:
     *  1. See if a binary archive is available, if so use that.
     *  2. See if an ascii archive is available, if so use that and make a binary one for next time.
     *  3. If nothing is available then the ascii archive is downloaded from
     *     http://www.cs.ox.ac.uk/people/gary.mirams/files/<model and pacing specific table>.tgz
     *     is downloaded, unpacked, loaded and converted to binary for next time.
     */
    void SetUpLookupTables();

public:
    /**
     * This constructor just sets some defaults.
     */
    ApPredictMethods();

    /**
     * Main running command.
     */
    virtual void Run();

    /**
     * Set the output directory
     *
     * @param rOuputDirectory  The directory to write results to - WILL BE WIPED!
     */
    void SetOutputDirectory(const std::string& rOuputDirectory);

    /**
     * This is a helper method to print out the available arguments.
     *
     * @return the arguments this class takes
     */
    static std::string PrintArguments();

    /**
     * @return The concentrations at which an action potential was evaluated.
     */
    std::vector<double> GetConcentrations(void);

    /**
     * @return The APD90s that were evaluated at the concentrations given by GetConcentrations().
     */
    std::vector<double> GetApd90s(void);

    /**
     * @return The 95% credible regions that are associated with the APD90 predictions given by
     * #GetApd90s().
     */
    std::vector<std::vector<double> > GetApd90CredibleRegions(void);
};

#endif //_APPREDICTMETHODS_HPP_
