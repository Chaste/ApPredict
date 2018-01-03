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

#include "AbstractDataStructure.hpp"
#include "Exception.hpp"

// Copied from https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
std::istream& AbstractDataStructure::SafeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for (;;)
    {
        int c = sb->sbumpc();
        switch (c)
        {
            case '\n':
                return is;
            case '\r':
                if (sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case EOF:
                // Also handle the case when the last line has no line ending
                if (t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}

void AbstractDataStructure::LoadDataFromFile(const std::string& rFileName, unsigned numHeaderLines)
{
    std::ifstream indata; // indata is like cin
    indata.open(rFileName.c_str()); // opens the file
    if (!indata.good())
    { // file couldn't be opened
        EXCEPTION("Couldn't open data file: " + rFileName);
    }

    bool first_line = true;
    unsigned num_lines_read = 0u;

    while (indata.good())
    {
        std::string this_line;
        SafeGetline(indata, this_line);
        num_lines_read++;

        if (this_line == "" || this_line == "\r")
        {
            if (indata.eof())
            { // If the blank line is the last line carry on OK.
                break;
            }
            else
            {
                EXCEPTION("No data found on line " << num_lines_read);
            }
        }
        std::stringstream line(this_line);

        if (first_line || (numHeaderLines > 0u && num_lines_read <= numHeaderLines))
        {
            first_line = false;
            // Try and read a header line if present
            if (LoadHeaderLine(line))
            {
                continue;
            }
        }

        // Load a standard data line.
        LoadALine(line);

        if (line.good())
        {
            std::string whats_left;
            line >> whats_left;
            EXCEPTION("These are unread items :'" << whats_left << "' on line " << num_lines_read << ", data reading structures may have bugs.");
        }
    }

    if (!indata.eof())
    {
        EXCEPTION("A file reading error occurred");
    }
}

bool AbstractDataStructure::LoadHeaderLine(std::stringstream& rLine)
{
    return false;
}
