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

#ifndef LOOKUPTABLELOADER_HPP_
#define LOOKUPTABLELOADER_HPP_

#include <boost/shared_ptr.hpp>

// ApPredict includes
#include "AbstractUntemplatedLookupTableGenerator.hpp"

/**
 * Class that will decide, based on the command line arguments supplied to ApPredict,
 * which is the most appropriate LookupTable to try and load. Most of the work is done
 * in the constructor which calls a bunch of private methods to decide which, if any,
 * lookup table is appropriate.
 *
 * The tables can be generated using TestMakeALookupTable.hpp
 *
 * Also looks at the manifest of ones that are available to download from the web and does
 * that if needed.
 */
class LookupTableLoader
{
private:
    /**
	   * The lookup table to return.
	   */
    boost::shared_ptr<AbstractUntemplatedLookupTableGenerator> mpLookupTable;

    /**
     * This is the set of channels that could be involved, and the ones that are,
     * in the order that they will appear in the filename of lookup tables
     * as generated by TestMakeALookupTable.hpp
     */
    std::vector<std::pair<std::string, bool> > mIdealChannelsInvolved;

    /** the CellML model's name (as specified by ODE GetSystemName() method). */
    std::string mModelName;

    /** the pacing frequency we are interested in */
    double mHertz;

    /**
     * The ideal lookup table.
     */
    std::string mIdealLookupTable;

    /**
     * The best available lookup table since mIdealLookupTable isn't available.
     */
    std::string mBestAvailableLookupTable;

    /**
     * Sets up #mIdealChannelsInvolved to match input arguments and sets #mIdealLookupTable string.
     */
    void DecideIdealTable();

    /**
     * A list of tables that could provide the information we need, in order of
     * preference (lowest dimension first).
     *
     * @return list of all the table base filenames that would work for this case.
     */
    std::vector<std::string> GenerateAllCompatibleTables();

    /**
     * A list of tables that are available at #mRemoteURL.
     *
     * Only returns tables that are relevant for this model and pacing rate.
     *
     * @return an unordered list of available tables.
     */
    std::vector<std::string> GetManifestOfTablesOnGarysWebsite();

    /**
     * A list of tables that are available in the current working directory.
     *
     * @return an unordered list of available tables.
     */
    std::vector<std::string> GetManifestOfLocalTablesInCwd();

    /**
     * Load this local file, convert to binary if we can
     */
    void LoadTableFromLocalBoostArchive(const std::string& rLookupTableBaseName);

    /**
     * Download and unzip a particular archive from #mRemoteURL.
     *
     * @param rArchiveFileBaseName  the archive to get.
     */
    void DownloadAndUnpack(const std::string& rArchiveFileBaseName);

    /**
     * The URL where the manifest and lookup tables are available from.
     * N.B. This string should be a folder and end in a "/".
     *
     * Value is set in the .cpp file.
     */
    static const std::string mRemoteURL;
     
public:
    /**
	   * Constructor
	   *
	   * @param rModelName  the CellML model's name (as specified by ODE GetSystemName() method).
	   * @param rHertz  the pacing frequency we are interested in.
	   */
    LookupTableLoader(const std::string& rModelName, const double& rHertz);

    /**
     * @return whether we have found and loaded a lookup table.
     */
    bool IsLookupTableAvailable();

    /**
     * @return the base filename for the ideal lookup table.
     */
    std::string GetIdealTable();

    /**
     * @return the base filename for the best lookup table if Ideal isn't available.
     */
    std::string GetBestAvailableTable();

    /**
	   * @return  a pointer to the lowest dimension LookupTable that is available.
	   */
    boost::shared_ptr<AbstractUntemplatedLookupTableGenerator> GetLookupTable();
};

#endif // LOOKUPTABLELOADER_HPP_
