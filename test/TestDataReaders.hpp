/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef TESTDATAREADERS_HPP_
#define TESTDATAREADERS_HPP_

#include <cxxtest/TestSuite.h>

#include "CardiovascRes2011DataStructure.hpp"

class TestDataReaders : public CxxTest::TestSuite
{
public:
    void TestDrugDataLoading(void) throw(Exception)
    {
        FileFinder file("projects/ApPredict/test/data/paper_drug_data.dat", RelativeTo::ChasteSourceRoot);

        // Test the drug data loads correctly...
        CardiovascRes2011DataStructure drug_data(file);
        TS_ASSERT_EQUALS(drug_data.GetNumDrugs(), 31u);

        unsigned ajmaline_idx = drug_data.GetDrugIndex("Ajmaline");
        TS_ASSERT_EQUALS(ajmaline_idx,0u);

        TS_ASSERT_THROWS_THIS(drug_data.GetDrugIndex("Sausages"),
                              "Drug Sausages not found.");

        unsigned quinidine = drug_data.GetDrugIndex("Quinidine");
        unsigned cisapride = drug_data.GetDrugIndex("Cisapride");
        unsigned tedisamil = drug_data.GetDrugIndex("Tedisamil");
        unsigned propranolol = drug_data.GetDrugIndex("Propranolol");
        unsigned verapamil = drug_data.GetDrugIndex("Verapamil");

        TS_ASSERT_EQUALS(drug_data.GetRedfernCategory(ajmaline_idx), 1u);
        TS_ASSERT_EQUALS(drug_data.GetDrugName(quinidine), "Quinidine");
        TS_ASSERT_DELTA(drug_data.GetIC50Value(quinidine,0), 16600, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetIC50Value(cisapride,2), 6.5, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetIC50Value(propranolol,1), 18000, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetIC50Value(propranolol,2), 2828, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetClinicalDoseRange(quinidine,0), 1000, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetClinicalDoseRange(quinidine,1), 4000, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetClinicalDoseRange(cisapride,0), 2, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetClinicalDoseRange(cisapride,1), 5, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetClinicalDoseRange(verapamil,0), 25, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetClinicalDoseRange(verapamil,1), 90, 1e-4);
        TS_ASSERT_EQUALS(drug_data.HasRedfernCategory(verapamil), true);
        TS_ASSERT_DELTA(drug_data.GetGrandiMeasure(verapamil), -20.753, 1e-4);

        TS_ASSERT_EQUALS(drug_data.HasClinicalDoseRange(verapamil), true);

        // We want a "no effect drug" to return DBL_MAX, which the CalculateConductanceFactor() method handles nicely (see below).
        unsigned chlorpromazine = drug_data.GetDrugIndex("Chlorpromazine");
        TS_ASSERT_EQUALS(drug_data.GetDrugName(chlorpromazine),"Chlorpromazine");
        TS_ASSERT_DELTA(drug_data.GetIC50Value(chlorpromazine,1), -2, 1e-9);

        // Check how it deals with a "NA" (No effect) entry - should return DBL_MAX for the IC50.
        // (i.e. a positive value which won't affect conductance so that analysis will run)
        TS_ASSERT_EQUALS(drug_data.GetDrugName(tedisamil), "Tedisamil");
        TS_ASSERT_DELTA(drug_data.GetIC50Value(tedisamil,0), 20000, 1e-9);
        TS_ASSERT_DELTA(drug_data.GetIC50Value(tedisamil,1), -2, 1e-9);
    }

    void TestConductanceFactorCalculations() throw(Exception)
    {
        // We've got two inputs that we want to return unchanged conductance:
        // -1 : we don't know what the conductance is
        // DBL_MAX : we know there is no effect.

        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(0.0,-1.0), 1.0, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(0.0,DBL_MAX), 1.0, 1e-9);

        // Test normal conductance calculations
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(1.0,1.0), 0.5, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(1.0,1.0,2.0), 0.5, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(2.0,1.0,2.0), 0.2, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(1.0,1.0,2.0,50), 0.75, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(0.0,1.0,2.0,90), 1.0, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(0.0,1.0,2.0,50), 1.0, 1e-9);

        // Saturation level (4th argument)
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(DBL_MAX,1.0,2.0,10), 0.1, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(DBL_MAX,1.0,2.0,50), 0.5, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(DBL_MAX,1.0,2.0,90), 0.9, 1e-9);

        TS_ASSERT_DELTA(AbstractDataStructure::ConvertIc50ToPic50(1000),3, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::ConvertPic50ToIc50(5),10, 1e-9);

        // Below here we are going to test some exceptional value handling
        // to deal with users doing strange things, for the web portal interface.
        TS_ASSERT_DELTA(AbstractDataStructure::ConvertPic50ToIc50( DBL_MAX), 0,       1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::ConvertPic50ToIc50(-DBL_MAX), DBL_MAX, 1e-9);

        TS_ASSERT_DELTA(AbstractDataStructure::ConvertIc50ToPic50(DBL_MAX), -302.2547, 1e-3);
        TS_ASSERT_DELTA(AbstractDataStructure::ConvertIc50ToPic50(0),        DBL_MAX,  1e-9);

        // Zero concentration implies unchanged conductance
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(0.0,0.0,1.0), 1.0, 1e-9);

        // Zero IC50 implies instantaneous block above here
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(0.001,0.0,1.0), 0.0, 1e-9);

        // Check that we get the also right response if we call CalculateConductanceFactor
        // with positive saturation - for agonists (activators) rather than antagonists (inhibitors).

        // Arguments are concentration, IC50, Hill, Saturation (defaults to 0 (%)).
        double saturation = 0.0;

        // Check usual inhibitor
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(1.0,1.0,1.0), 0.5, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(1.0,1.0,1.0,saturation), 0.5, 1e-9);

        // Now do an agonist
        saturation = 150.0;

        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(1.0,1.0,1.0,saturation), 1.25, 1e-9);

        saturation = 125.0;
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(1.0,1.0,1.0,saturation), 1.125, 1e-9);
        TS_ASSERT_DELTA(AbstractDataStructure::CalculateConductanceFactor(DBL_MAX,1.0,1.0,saturation), 1.25, 1e-9);
    }
};

#endif // TESTDATAREADERS_HPP_
