/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TORSADEPREDICTMETHODS_HPP_
#define TORSADEPREDICTMETHODS_HPP_

#include "ApPredictMethods.hpp"
#include "LinearDiscriminantAnalysis.hpp"

/**
 * This class calls ApPredict to predict action potentials
 * for given ion channel IC50 and Hill data. For a range of
 * concentrations as requested. It then performs a linear
 * discrimanant analysis to predict the Torsade risk category
 * of your drug at each concentration, as per Mirams et al. 2011.
 *
 * Call PrintArguments to see usage.
 */
class TorsadePredictMethods : public ApPredictMethods
{
public:
    /**
     * This constructor does all the work of running the ActionPotential (AP) preDiCT
     * with Torsadogenic risk estimates (TorsadePredict) program
     */
    TorsadePredictMethods();

    /**
     * Run the predictions.
     */
    virtual void Run();

    /**
     * This is a helper method to print out the available arguments.
     *
     * @return available arguments for this class.
     */
    static std::string PrintArguments();

    /**
     * Get the predicted Redfern category for each concentration of this compound.
     *
     * @return The torsade predictions.
     */
    std::vector<unsigned> GetTorsadePredictions();

protected:
    /**
     * Load up the Torsade training data from file.
     * This relies on "paper_drug_data.dat" being present.
     *
     * @return a Linear Discriminant Analysis (LDA) object preloaded with the Grandi APD90 training data.
     */
    LinearDiscriminantAnalysis LoadLdaFromDrugData();

    /**
     * Use the APD data from the current run to
     * call the LDA methods and classify points.
     */
    void MakeTorsadePredictions();

private:
    friend class TestTorsadePredict;

    /** A vector used to store the Torsade predictions (as indexed by Redfern i.e. 0,1,2,3 = cats 2/3/4/5) */
    std::vector<unsigned> mTorsadePredictions;

    /**
     * Write the TdP classification results to HTML table.
     */
    void WriteTorsadeResultsToFile();

};

#endif // TORSADEPREDICTMETHODS_HPP_
