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

#include "TorsadePredictMethods.hpp"
#include "SetupModel.hpp"
#include "CardiovascRes2011DataStructure.hpp"
#include "FileFinder.hpp"

std::string TorsadePredictMethods::PrintArguments()
{
    std::string message = "TorsadePredict::Please provide these inputs:\n";
    message += PrintCommonArguments();
    return message;
}

TorsadePredictMethods::TorsadePredictMethods()
 : ApPredictMethods()
{
}

void TorsadePredictMethods::Run()
{
    mProgramName = "Torsade PreDiCT";
    mOutputFolder = "TorsadePredict_output";

    SetupModel setup(this->mHertz, 5u); // Hardcoded to Grandi model.
    mpModel = setup.GetModel();

    CommonRunMethod();

    assert(mComplete);

    // After calling the normal methods
    // Additional steps for Torsade runs included here.
    MakeTorsadePredictions();

    // This is in a separate method to the above so we can test the above on its own...
    WriteTorsadeResultsToFile();
}

void TorsadePredictMethods::MakeTorsadePredictions()
{
    if (mApd90s.size()==0)
    {
        EXCEPTION("APDs do not appear to have been recorded.");
    }
    mTorsadePredictions.reserve(mApd90s.size());

    // Work out the Measure
    std::vector<double> largest_percent_change;
    largest_percent_change.push_back(0.0); // First is always control
    for (unsigned i=1; i<mApd90s.size(); i++)
    {
        double percent = 100*(mApd90s[i]-mApd90s[0])/mApd90s[0];
        largest_percent_change.push_back(percent);
    }

    // Second pass to check for largest +ve effect at low dose...
    // This makes the measure "Largest Effect Dose" rather than just "dose".
    for (unsigned i=1; i<mApd90s.size(); i++)
    {
        // Check over lower doses
        for (unsigned j=0; j<i; j++)
        {
            // If the lower dose made the AP longer than control (and longer than this dose) then replace this APD with that one.
            if (largest_percent_change[j] > 0 && largest_percent_change[j] > largest_percent_change[i])
            {
                largest_percent_change[i] = largest_percent_change[j];
            }
        }
    }

    // Perform the LDA
    LinearDiscriminantAnalysis lda = LoadLdaFromDrugData();
    for (unsigned i=0; i<mApd90s.size(); i++)
    {
        vector<double> test_point(1);
        test_point(0) = largest_percent_change[i];
        mTorsadePredictions.push_back(lda.ClassifyThisPoint(test_point)+2u); // We add two because our redfern categories start at 2 not 0.
    }
}

void TorsadePredictMethods::WriteTorsadeResultsToFile()
{
    assert(mpFileHandler);

    // Open an output file for the Torsade results
    out_stream torsade_results_file = mpFileHandler->OpenOutputFile("tdp_results.html");
    *torsade_results_file << "<html>\n<head><title>Torsade preDiCT Results</title></head>\n";
    *torsade_results_file << "<STYLE TYPE=\"text/css\">\n<!--\nTD{font-size: 12px;}\n--->\n</STYLE>\n";
    *torsade_results_file << "<body>\n";
    *torsade_results_file << "<table width=\"60%\" style=\"background-color:white\" border=\"1\" cellpadding=\"2\" cellspacing=\"0\">\n";
    *torsade_results_file << "<tr><td>Concentration (uM)</td><td>APD90 (ms)</td><td>Risk Category Prediction</td></tr>\n"; // Header line
    std::vector<std::string> colours;
    colours.push_back("Red");
    colours.push_back("Orange");
    colours.push_back("Limegreen");
    colours.push_back("Limegreen");
    for (unsigned i=0; i<mTorsadePredictions.size(); i++)
    {
        if (!mSuppressOutput) std::cout << "Conc = " << mConcs[i] << "uM, APD90 = " << mApd90s[i] << "ms, risk prediction =  " << mTorsadePredictions[i] << "\n";// << std::flush;
        *torsade_results_file << "<tr style=\"background-color:" << colours[mTorsadePredictions[i]-2] << "\"><td>"<< mConcs[i] << "</td><td>" << mApd90s[i] << "</td><td>" << mTorsadePredictions[i] << "</td></tr>\n";
    }
    *torsade_results_file << "</table>\n</body>\n</html>\n";
    torsade_results_file->close();
}


LinearDiscriminantAnalysis TorsadePredictMethods::LoadLdaFromDrugData()
{
    // We have to look for the drug data in a couple of places because the executable won't have the Chaste source.
    const std::string drug_data_file = "paper_drug_data.dat";

    boost::shared_ptr<FileFinder> p_file_finder = boost::shared_ptr<FileFinder>(new FileFinder("projects/ApPredict/test/data/" + drug_data_file, RelativeTo::ChasteSourceRoot));
    if (!p_file_finder->Exists())
    {
        p_file_finder->SetPath(drug_data_file, RelativeTo::CWD);
        if (!p_file_finder->Exists())
        {
            EXCEPTION("The file \"" << drug_data_file << "\" should be in the current working directory and is missing.");
        }
    }
    CardiovascRes2011DataStructure drug_data(*p_file_finder);

    matrix<double> cat1and2 (1,1);
    matrix<double> cat3 (1,1);
    matrix<double> cat4 (1,1);
    matrix<double> cat5 (1,1);
    cat1and2(0,0) = cat3(0,0) = cat4(0,0) = cat5(0,0) = DBL_MAX;

    // Generate structures for training data...
    for (unsigned i=0; i<drug_data.GetNumDrugs(); i++)
    {
        double grandi_measure;
        try
        {
            grandi_measure = drug_data.GetGrandiMeasure(i);
        }
        catch (Exception& e)
        {
            assert(e.GetShortMessage()=="No data available on Grandi measure for this drug.");
            continue;
        }

        if (drug_data.GetDrugName(i)=="Ranolazine")
        {   // We aren't currently using Ranolazine for the training.
            continue;
        }

        unsigned redfern = drug_data.GetRedfernCategory(i);
        if (redfern==1 || redfern==2)
        {
            if (cat1and2(0,0) != DBL_MAX)
            {
                cat1and2.resize(cat1and2.size1()+1,1);
            }
            cat1and2(cat1and2.size1()-1,0) = grandi_measure;
        }
        else if (redfern==3)
        {
            if (cat3(0,0) != DBL_MAX)
            {
                cat3.resize(cat3.size1()+1,1);
            }
            cat3(cat3.size1()-1,0) = grandi_measure;
        }
        else if (redfern==4)
        {
            if (cat4(0,0) != DBL_MAX)
            {
                cat4.resize(cat4.size1()+1,1);
            }
            cat4(cat4.size1()-1,0) = grandi_measure;
        }
        else
        {
            assert(redfern==5);
            if (cat5(0,0) != DBL_MAX)
            {
                cat5.resize(cat5.size1()+1,1);
            }
            cat5(cat5.size1()-1,0) = grandi_measure;
        }
    }

    std::vector<matrix<double> > training;
    training.push_back(cat1and2);
    training.push_back(cat3);
    training.push_back(cat4);
    training.push_back(cat5);

    LinearDiscriminantAnalysis lda(training);
    return lda;
}


std::vector<unsigned> TorsadePredictMethods::GetTorsadePredictions(void)
{
    if (!mComplete)
    {
        EXCEPTION("Simulation has not been run - check arguments.");
    }
    return mTorsadePredictions;
}


