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

#ifndef ABSTRACTDATASTRUCTURE_HPP_
#define ABSTRACTDATASTRUCTURE_HPP_

#include <fstream>
#include <vector>

#include "UblasVectorInclude.hpp"
#include "FileFinder.hpp"

/**
 * A class which is designed to read in data from a text file, subclassed for particular formats.
 */
class AbstractDataStructure
{
protected:

    /**
     * Main method for loading the data file, line by line it calls
     * the overloaded method LoadALine
     *
     * @param fileName  of the drug data file
     */
    void LoadDataFromFile(std::string fileName);

    /**
     * This method must be overridden by subclasses to fill in the
     * data structures with required info from file
     *
     * @param rLine  a line in stringsteam format.
     */
    virtual void LoadALine(std::stringstream& rLine)=0;

    /**
     * Read a header line if present.
     *
     * If not implemented this class doesn't touch the line, and returns false.
     * Subclasses should implement the line loading and return true.
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

};

#endif // ABSTRACTDATASTRUCTURE_HPP_
