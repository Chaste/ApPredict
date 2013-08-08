/*

Copyright (c) 2005-2013, University of Oxford.
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
#include "OutputFileHandler.hpp"

/**
 * Common code to allow this to be run as a test via scons and also as a
 * executable application from the apps folder.
 *
 * The member variables are only used by the Torsade prediction methods.
 */
class ApPredictMethods : public AbstractActionPotentialMethod
{
private:
    /**
     * A method to avoid lots of copying and pasting in the main drug block application method.
     *
     * @param pModel  The CVODE model.
     * @param rMetadataName  The Oxford metadata name for the conductance we want to block.
     * @param rShortName  The short version of the name for screen output.
     * @param default_conductance  The default value of the conductance before we started messing.
     * @param concentration  The current drug concentration.
     * @param iC50  The IC50 for this channel.
     * @param hill  The Hill for this channel.
     */
    void ApplyDrugBlock(boost::shared_ptr<AbstractCvodeCell> pModel,
                        const std::string& rMetadataName,
                        const std::string& rShortName,
                        const double default_conductance,
                        const double concentration,
                        const double iC50,
                        const double hill);

protected:

    /** Whether the simulation completed successfully */
    bool mComplete;

    /** A vector used to store the APD90s calculated in the main method */
    std::vector<double> mAPD90s;

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
     * @param rIc50  a default IC50 value (usually -1), overwritten if argument is present to a value (assumed to be in uM).
     * @param rHill  a default Hill coefficient (usually -1), overwritten if argument is present.
     * @param rChannel  The name of the channel in the input arguments
     */
    void ReadInIC50AndHill(double& rIc50, double& rHill, const std::string& rChannel);

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
     * @return the arguments both ApPredict and TorsadePredict take.
     */
    static std::string PrintCommonArguments();

    /**
     * A run method common to both ApPredict and TorsadePredict.
     *
     * Does the majority of the work!
     */
    void CommonRunMethod();

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

};

#endif //_APPREDICTMETHODS_HPP_
