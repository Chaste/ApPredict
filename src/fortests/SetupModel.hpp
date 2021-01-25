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

#ifndef SETUPMODEL_HPP_
#define SETUPMODEL_HPP_

#include <boost/shared_ptr.hpp>

#include "AbstractCvodeCell.hpp"
#include "OutputFileHandler.hpp"

/**
 * Class to return a Cvode cell model with appropriate stimulus based on the
 * CellML default stimulus. This class reads the command line argument
 * --model <unsigned> according to the options in constructor documentation
 * if no model is specifically requested.
 */
class SetupModel
{
private:
    /** Private default constructor to stop this being called and point in the direction of the other constructor */
    SetupModel(){};

    /** OutputFileHander to specify a directory to use when autogenerating CellML */
    boost::shared_ptr<OutputFileHandler> mpHandler;

    /** The cell model (in CVODE format) */
    boost::shared_ptr<AbstractCvodeCell> mpModel;

public:
    /**
     * This method generates a cell model for one of the models specified by #model_index.
     * If model_index is not set then this class attempts to read it from CommandLineOptions as
     * specified in #PrintArguments().
     *
     * It attempts to use the model's default stimulus magnitude and duration from CellML
     * with a frequency of #hertz and a start time of 1ms.
     *
     * If this is not present it uses sensible defaults.
     *
     * @param model_index  1 = Shannon, 2=TenTusscher, 3 = Mahajan, 4 = Hund-Rudy, 5 = Grandi,
     *        6 = O'Hara-Rudy 2011, 7 = Paci ventricular, 8 = O'Hara-Rudy CiPA v1 2017, 9 = Faber-Rudy
     * @param hertz  The frequency of the regular stimulus that this model should use.
     * @param pHandler  An optional pointer to use as a working directory when
     *                  generating code on the fly from CellML (defaults to empty pointer).
     */
    SetupModel(const double& rHertz,
               unsigned model_index = UNSIGNED_UNSET,
               boost::shared_ptr<OutputFileHandler> pHandler = boost::shared_ptr<OutputFileHandler>());

    /**
     * Print the option and possible arguments that this class takes.
     *
     * @return arguments for command line display.
     */
    static std::string PrintArguments()
    {
        return "* EITHER --model\n"
               "*   options: 1 = Shannon, 2 = TenTusscher (06), 3 = Mahajan,\n"
               "*            4 = Hund-Rudy, 5 = Grandi, 6 = O'Hara-Rudy 2011 (endo),\n"
               "*            7 = Paci (ventricular), 8 = O'Hara-Rudy CiPA v1 2017 (endo)\n"
               "*            9 = Faber-Rudy.\n"
               "* OR --cellml <file>\n";
    }

    /**
     * Get a boost shared pointer to the model the constructor sets up
     * @return a pointer to the model
     */
    boost::shared_ptr<AbstractCvodeCell> GetModel();
};

#endif // SETUPMODEL_HPP_
