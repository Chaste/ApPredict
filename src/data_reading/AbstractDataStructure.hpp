/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef ABSTRACTDATASTRUCTURE_HPP_
#define ABSTRACTDATASTRUCTURE_HPP_

#include <fstream>
#include <vector>

#include "Exception.hpp"
#include "FileFinder.hpp"
#include "UblasVectorInclude.hpp"

/**
 * A class which is designed to read in data from a text file, subclassed for particular formats.
 *
 * Also provides some helper methods for calculating conductance block and converting from pIC50 to IC50s.
 */
class AbstractDataStructure
{
private:
    std::istream& SafeGetline(std::istream& is, std::string& t);

protected:
    /**
     * Main method for loading the data file, line by line it calls
     * the overloaded method LoadALine
     *
     * @param rFileName  of the drug data file
     * @param numHeaderLines  the number of header lines in the file to skip (your class MUST implement LoadHeaderLine to use this although it can just return true).
     */
    void LoadDataFromFile(const std::string& rFileName, unsigned numHeaderLines = 0);

    /**
     * This method must be overridden by subclasses to fill in the
     * data structures with required info from file
     *
     * @param rLine  a line in stringsteam format.
     */
    virtual void LoadALine(std::stringstream& rLine) = 0;

    /**
     * Read a header line if present.
     *
     * If not implemented this class doesn't touch the line, and returns false (no header is considered present)
     * Subclasses with a header line should implement the line loading and return true.
     */
    virtual bool LoadHeaderLine(std::stringstream& rLine);

public:
    /**
     * Default Constructor (empty)
     */
    AbstractDataStructure(){};

    /**
     * Destructor (empty)
     */
    virtual ~AbstractDataStructure(){};

    /**
     * Calculate the probability of a channel being open given this drug, IC50 and hill coefficient.
     *
     * Note: A negative IC50 value is interpreted as "drug has no effect on this channel".
     *
     * @param rConc  concentration of the drug.
     * @param rIC50  IC50 value for this drug and channel (note if saturation is less than 100 this parameter is not IC50, but 50% of max  response)
     * @param hill  Hill coefficient for this drug dependent inactivation curve (defaults to 1).
     * @param saturation  The saturation level for this drug (defaults to 0), relative to unchanged channel conductance at 100%.
     *
     * @return proportion of channels which are still active at this drug concentration
     */
    static double CalculateConductanceFactor(const double& rConc,
                                             const double& rIC50,
                                             double hill = 1.0,
                                             double saturation = 0)
    {
        if (saturation < 0)
        {
            // default is a full inhibitor - saturates at zero conductance.
            saturation = 0.0;
        }

        // To avoid divide-by-zero style stuff, if there is no drug, conductance must be unchanged.
        if (rConc == 0)
        {
            return 1.0;
        }

        if (rIC50 < 0) // missing information ('-1'), or known-to-be-no-effect ('-2'), do not alter conductance.
        {
            return 1.0;
        }

        // If the hill coefficient has not been set (defaults to a negative value) then use 1.0 instead.
        if (hill < 0)
        {
            hill = 1.0;
        }

        return 1.0 - ((100.0 - saturation) / 100.0) * (1.0 - 1.0 / (1.0 + pow((rConc / rIC50), hill)));
    }

    /**
     * Converts an IC50 IN MICRO MOLAR (uM) into a pIC50 (in log Molar)
     * @param rIc50  The IC50 value in microMolar
     * @return the pIC50 value
     */
    static double ConvertIc50ToPic50(const double& rIc50)
    {
        double result = -log10((1e-6) * rIc50);
        // Handle very small IC50 values gracefully (result = -Inf)
        if (!std::isfinite(result))
        {
            result = DBL_MAX;
        }
        return result;
    }

    /**
     * Converts a pIC50 value into an IC50 IN MICRO MOLAR (uM)
     * @return rIc50  The IC50 value in microMolar
     * @param the pIC50 value
     */
    static double ConvertPic50ToIc50(const double& rPic50)
    {
        double result = pow(10.0, 6.0 - rPic50);

        // Handle large negative IC50 values gracefully (result = -Inf or Inf)
        if (!std::isfinite(result))
        {
            result = DBL_MAX;
        }
        return result;
    }
};

#endif // ABSTRACTDATASTRUCTURE_HPP_
