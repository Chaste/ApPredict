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

#ifndef _TESTDOSECALCULATOR_HPP_
#define _TESTDOSECALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "DoseCalculator.hpp"


class TestDoseCalculator : public CxxTest::TestSuite
{
public:
    void TestDoseCalculatorSpecifyConcs(void) throw(Exception)
    {
		// This tests the simpler interface
		std::vector<double> concs;
		concs.push_back(1.1);
		concs.push_back(1.5);
		concs.push_back(1.6);
		concs.push_back(1.7);

		DoseCalculator calculator(concs);

		std::vector<double> test_concs = calculator.GetConcentrations();

		// We should have one more test concentration (control added)
		TS_ASSERT_EQUALS(test_concs.size(), concs.size() + 1u);

		TS_ASSERT_DELTA(test_concs[0], 0.0, 1e-9);

		for (unsigned i=0; i<concs.size() ; i++)
		{
		 TS_ASSERT_DELTA(test_concs[i+1], concs[i], 1e-9);
		}

		// Now setup for a log scale (should add 1nM as well as 0nM)
		calculator.SetLogScale(true);
		test_concs = calculator.GetConcentrations();

		TS_ASSERT_EQUALS(test_concs.size(), concs.size() + 2u);

		TS_ASSERT_DELTA(test_concs[0], 0.0, 1e-9);
		TS_ASSERT_DELTA(test_concs[1], 1e-3, 1e-9);

		for (unsigned i=0; i<concs.size() ; i++)
		{
			TS_ASSERT_DELTA(test_concs[i+2], concs[i], 1e-9);
		}
	}

	void TestDoseCalculatorForInsertedSubValues(void) throw(Exception)
    {
		{ // For non-log scale
			std::vector<double> concs;
			concs.push_back(0);
			concs.push_back(10);
			concs.push_back(20);
			concs.push_back(30);

			DoseCalculator calculator(concs);
			calculator.SetNumSubdivisions(1);

			std::vector<double> result_concs = calculator.GetConcentrations();

			TS_ASSERT_EQUALS(result_concs.size(), 7u);
			TS_ASSERT_DELTA(result_concs[0],0,1e-9);
			TS_ASSERT_DELTA(result_concs[1],5,1e-9);
			TS_ASSERT_DELTA(result_concs[2],10,1e-9);
			TS_ASSERT_DELTA(result_concs[3],15,1e-9);
			TS_ASSERT_DELTA(result_concs[4],20,1e-9);
			TS_ASSERT_DELTA(result_concs[5],25,1e-9);
			TS_ASSERT_DELTA(result_concs[6],30,1e-9);
		}

		{ // For log scale
			std::vector<double> concs;
			concs.push_back(0);
			concs.push_back(0.001);
			concs.push_back(0.01);
			concs.push_back(0.1);

			DoseCalculator calculator(concs);
			calculator.SetNumSubdivisions(1);
			calculator.SetLogScale(true);

			std::vector<double> result_concs = calculator.GetConcentrations();

			TS_ASSERT_EQUALS(result_concs.size(), 6u);
			TS_ASSERT_DELTA(result_concs[0],0		    ,1e-9);
			TS_ASSERT_DELTA(result_concs[1],0.001       ,1e-9);
			TS_ASSERT_DELTA(result_concs[2],0.00316228	,1e-7);
			TS_ASSERT_DELTA(result_concs[3],0.01		,1e-9);
			TS_ASSERT_DELTA(result_concs[4],0.0316228	,1e-7);
			TS_ASSERT_DELTA(result_concs[5],0.1		    ,1e-9);
		}

		{ // For log scale
			std::vector<double> concs;
			concs.push_back(0);
			concs.push_back(1);
			concs.push_back(10);
			concs.push_back(100);

			DoseCalculator calculator(concs);
			calculator.SetNumSubdivisions(1);
			calculator.SetLogScale(true);

			std::vector<double> result_concs = calculator.GetConcentrations();

			TS_ASSERT_EQUALS(result_concs.size(), 8u);
			TS_ASSERT_DELTA(result_concs[0],0		,1e-9);
			TS_ASSERT_DELTA(result_concs[1],1e-3	,1e-9);
			TS_ASSERT_DELTA(result_concs[2],0.0316228,1e-6);
			TS_ASSERT_DELTA(result_concs[3],1       ,1e-9);
			TS_ASSERT_DELTA(result_concs[4],3.16228	,1e-5);
			TS_ASSERT_DELTA(result_concs[5],10		,1e-9);
			TS_ASSERT_DELTA(result_concs[6],31.6228	,1e-4);
			TS_ASSERT_DELTA(result_concs[7],100		,1e-9);
		}

		{ // For log scale
			std::vector<double> concs;
			concs.push_back(30);
			concs.push_back(100);
			concs.push_back(300);
			concs.push_back(1000);

			DoseCalculator calculator(concs);
			calculator.SetNumSubdivisions(2);
			calculator.SetLogScale(true);

			std::vector<double> result_concs = calculator.GetConcentrations();

			TS_ASSERT_EQUALS(result_concs.size(), 14u);
			TS_ASSERT_DELTA(result_concs[0],0		,1e-9);
			TS_ASSERT_DELTA(result_concs[1],1e-3	,1e-9);
			TS_ASSERT_DELTA(result_concs[4],30		,1e-9);
			TS_ASSERT_DELTA(result_concs[7],100		,1e-9);
			TS_ASSERT_DELTA(result_concs[10],300	,1e-9);
			TS_ASSERT_DELTA(result_concs[13],1000	,1e-9);
		}
	}

	void TestForSpreadOfValues(void) throw(Exception)
	{
		TS_ASSERT_THROWS_THIS(DoseCalculator calc(0,1),
				"Bottom test concentration cannot be larger than top test concentration.");

		// For equally spaced
		{
			DoseCalculator calc(1000,0);
			std::vector<double> result_concs = calc.GetConcentrations();

			TS_ASSERT_EQUALS(result_concs.size(), 11u);
			for (unsigned i=0; i<result_concs.size(); i++)
			{
				TS_ASSERT_DELTA(result_concs[i], ((double)(i))*100, 1e-9);
			}
		}
	}

    void TestForSpreadStartingAboveZero(void) throw(Exception)
    {
        // For equally spaced starting above zero
        DoseCalculator calc(1000,500);
        calc.SetNumSubdivisions(9u);
        std::vector<double> result_concs = calc.GetConcentrations();

        TS_ASSERT_EQUALS(result_concs.size(), 12u);

        // Check it adds in a control
        TS_ASSERT_DELTA(result_concs[0], 0.0, 1e-9);
        for (unsigned i=1; i<result_concs.size(); i++)
        {
            TS_ASSERT_DELTA(result_concs[i], 500+(((double)(i))-1)*50, 1e-9);
        }
    }


	// For log-spaced starting at zero
    void TestForLogSpacedStartingAtZero(void) throw(Exception)
    {
        DoseCalculator calc(1000,0);
        calc.SetLogScale(true);
        calc.SetNumSubdivisions(5u);
        std::vector<double> result_concs = calc.GetConcentrations();

        // 8 = 5 subdivisions + top + bottom + one
        TS_ASSERT_EQUALS(result_concs.size(), 8u);
        TS_ASSERT_DELTA(result_concs[0], 0, 1e-9);
        TS_ASSERT_DELTA(result_concs[1], 1e-3, 1e-9);
        for (unsigned i=2; i<result_concs.size(); i++)
        {
            TS_ASSERT_DELTA(result_concs[i], pow(10,(((double)(i))-4)), 1e-9);
        }
    }

    // For log-spaced starting at 1nM
    void TestForLogSpacedStartingAtOneNanoMolar(void) throw(Exception)
    {
        DoseCalculator calc(1,1e-3);
        calc.SetLogScale(true);
        calc.SetNumSubdivisions(5u);
        std::vector<double> result_concs = calc.GetConcentrations();

        // 8 = 5 subdivisions + top + bottom + one
        TS_ASSERT_EQUALS(result_concs.size(), 8u);
        TS_ASSERT_DELTA(result_concs[0], 0, 1e-9);
        TS_ASSERT_DELTA(result_concs[1], 1e-3, 1e-9);
        for (unsigned i=2; i<result_concs.size(); i++)
        {
            TS_ASSERT_DELTA(result_concs[i], pow(10,((((double)(i))-1)/2) - 3), 1e-9);
        }
    }

    // For log-spaced starting above 1nM
    void TestForLogSpacedStartingAboveOneNanoMolar(void) throw(Exception)
    {
        DoseCalculator calc(1000,100);
        calc.SetLogScale(true);
        calc.SetNumSubdivisions(1u);
        std::vector<double> result_concs = calc.GetConcentrations();

        // 5 = 2 subdivisions (one in range, one between bottom and one) + top + bottom + one nM + control
        TS_ASSERT_EQUALS(result_concs.size(), 5u);
        TS_ASSERT_DELTA(result_concs[0], 0,        1e-9);
        TS_ASSERT_DELTA(result_concs[1], 1e-3,     1e-9);
        TS_ASSERT_DELTA(result_concs[2], 100,      1e-9);
        TS_ASSERT_DELTA(result_concs[3], 316.228,  1e-3);
        TS_ASSERT_DELTA(result_concs[4], 1000,     1e-9);
    }

    // For log-spaced starting above 1nM with extra spacing in first position
    void TestForLogSpacedStartingAtOneMicroMolarExtraSpacing(void) throw(Exception)
    {
        DoseCalculator calc(10000,1);
        calc.SetLogScale(true);
        calc.SetNumSubdivisions(1u);
        std::vector<double> result_concs = calc.GetConcentrations();

        // 5 = 1 subdivisions (one in range) + top + bottom + one + control
        TS_ASSERT_EQUALS(result_concs.size(), 5u);
        TS_ASSERT_DELTA(result_concs[0], 0,         1e-9);
        TS_ASSERT_DELTA(result_concs[1], 1e-3,      1e-9);
        TS_ASSERT_DELTA(result_concs[2], 1,         1e-9);
        TS_ASSERT_DELTA(result_concs[3], 100,       1e-4);
        TS_ASSERT_DELTA(result_concs[4], 10000,     1e-9);
    }

    // For log-spaced starting above 1nM with extra spacing in first position
    void TestForLogSpacedStartingAboveOneNanoMolarExtraSpacing(void) throw(Exception)
    {
        DoseCalculator calc(10000,1000);
        calc.SetLogScale(true);
        calc.SetNumSubdivisions(3u);
        std::vector<double> result_concs = calc.GetConcentrations();

        // 10 = 3x2 subdivisions (one in range, one between bottom and one) + top + bottom + one + control
        TS_ASSERT_EQUALS(result_concs.size(), 7u);
        TS_ASSERT_DELTA(result_concs[0], 0,         1e-9);
        TS_ASSERT_DELTA(result_concs[1], 1e-3,      1e-9);
        TS_ASSERT_DELTA(result_concs[2], 1000,      1e-9);
        TS_ASSERT_DELTA(result_concs[3], 1778.28,   1e-2);
        TS_ASSERT_DELTA(result_concs[4], 3162.28,   1e-2);
        TS_ASSERT_DELTA(result_concs[5], 5623.41,   1e-2);
        TS_ASSERT_DELTA(result_concs[6], 10000,     1e-9);
    }

    // For log-spaced starting above 1nM with extra spacing in first position
    void TestForLogSpaced100NanoMolar(void) throw(Exception)
    {
        DoseCalculator calc(10,0.1);
        calc.SetLogScale(true);
        calc.SetNumSubdivisions(3u);
        std::vector<double> result_concs = calc.GetConcentrations();

        // 10 = 3x2 subdivisions (one in range, one between bottom and one) + top + bottom + one + control
        TS_ASSERT_EQUALS(result_concs.size(), 7u);
        TS_ASSERT_DELTA(result_concs[0], 0,         1e-9);
        TS_ASSERT_DELTA(result_concs[1], 1e-3,      1e-9);
        TS_ASSERT_DELTA(result_concs[2], 0.1,       1e-7);
        TS_ASSERT_DELTA(result_concs[3], 0.316228,  1e-6);
        TS_ASSERT_DELTA(result_concs[4], 1,         1e-4);
        TS_ASSERT_DELTA(result_concs[5], 3.16228,   1e-2);
        TS_ASSERT_DELTA(result_concs[6], 10,        1e-9);
    }

    // For a corner-case Geoff found where bottom and top doses are
    // set less than control on log scale (1nM).
    void TestForLogSpacedVeryLowRange(void) throw(Exception)
    {
        // 0.0001 uM == 0.1nM which is < control of 1nM for log scale.
        DoseCalculator calc(0.0002,0.0001);
        calc.SetLogScale(true);
        calc.SetNumSubdivisions(1u);
        std::vector<double> result_concs = calc.GetConcentrations();

        TS_ASSERT_EQUALS(result_concs.size(), 4u);
        TS_ASSERT_DELTA(result_concs[0], 0,           1e-9);
        TS_ASSERT_DELTA(result_concs[1], 0.0001,      1e-9);
        TS_ASSERT_DELTA(result_concs[2], 0.000141421, 1e-9);
        TS_ASSERT_DELTA(result_concs[3], 0.0002,      1e-6);
    }

};

#endif // _TESTDOSECALCULATOR_HPP_

