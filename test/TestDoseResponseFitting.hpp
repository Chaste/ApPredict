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

#ifndef TESTDOSERESPONSEFITTING_HPP_
#define TESTDOSERESPONSEFITTING_HPP_

#include <boost/assign/list_of.hpp>
#include <cassert>
#include <cmath>
#include <cxxtest/TestSuite.h>
#include <iostream>
#include <string>
#include <vector>

#include "HillFunction.hpp"
#include "NelderMeadMinimizer.hpp"
#include "OutputFileHandler.hpp"
#include "RunHillFunctionMinimization.hpp"
#include "Timer.hpp"
#include "UblasVectorInclude.hpp"

// Main function running code for test compound data to check that known
// parameters are returned
class TestDoseResponseFitting : public CxxTest::TestSuite
{
private:
    void PrintParameterValues(std::vector<double> parameters)
    {
        std::cout << "Parameters size = " << parameters.size() << "\n";
        std::cout << "Fitted parameters: \n";
        for (unsigned i = 0; i < parameters.size(); i++)
        {
            std::cout << "Parameter[" << i << "] = " << parameters[i] << "\n";
        }
    }

    // Put in some timings to work out the tricky cases...
    void setUp() { Timer::Reset(); }
    void tearDown() { Timer::Print("Test Elapsed"); }

    std::vector<double> CallDoseResponseFittingCode(
        const std::vector<double>& rDoses, const std::vector<double>& rResponses,
        unsigned numToFit = 2u, bool roundValues = false)
    {
        // Check that at least 1 data point and the same number of concentration and
        // inhibition data points are input
        TS_ASSERT(rDoses.size() == rResponses.size());
        TS_ASSERT(rDoses.size() > 0);

        // Runs minimization algorithm for given input data
        RunHillFunctionMinimization testcompound(rDoses, rResponses, numToFit,
                                                 roundValues);
        return testcompound
            .Run(); // this returns the parameters that have been fitted.
    }

public:
    void TestDoseResponseFittingPoints1()
    {
        // Concentration and inhibition data for compound input for test case 1 -
        // fitting with
        // just one data point input
        std::vector<double> testconcentrations = boost::assign::list_of(10);
        std::vector<double> testinhibitions = boost::assign::list_of(50);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 1u);

        // Prints to screen whether first test case has been passed - checks that
        // correct
        // parameter is returned and that only one parameter is returned

        if (parameters[0] == 10 && parameters.size() == 1)
        {
            std::cout << "Test Case 1 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 1 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints2()
    {
        // Concentration and inhibition data for compound input for test case 2 -
        // fitting with
        // more than one data point input
        std::vector<double> testconcentrations = boost::assign::list_of(0.37)(1.11)(3.33)(10);
        std::vector<double> testinhibitions = boost::assign::list_of(27.49)(51.45)(74.8)(88.49);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        // Checks that know parameter results are returned and that exactly 2
        // parameters are returned
        // and prints whether test case has been passed to screen
        if ((fabs(parameters[0] - 1.046) < 1e-3) && (fabs(parameters[1] - 0.9250) < 1e-3) && parameters.size() == 2)
        {
            std::cout << "Test Case 2 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 2 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }

        parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 1u);

        // Checks that know parameter results are returned and that exactly 1
        // parameters are returned
        // and prints whether test case has been passed to screen
        if ((fabs(parameters[0] - 1.0581) < 1e-4) && parameters.size() == 1u)
        {
            std::cout << "Test Case 2a Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 2a Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints3()
    {
        // Concentration and inhibition data for compound input for test case 3 -
        // fitting with
        // more than one data point input
        std::vector<double> testconcentrations = boost::assign::list_of(0.37)(1.11)(3.33)(10);
        std::vector<double> testinhibitions = boost::assign::list_of(7.0727)(17.61178)(37.5152)(62.7956);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        // Checks that know parameter results are returned and that exactly 2
        // parameters are returned
        // and prints whether test case has been passed to screen
        if ((fabs(parameters[0] - 5.729) < 1e-3) && (fabs(parameters[1] - 0.940) < 1e-3) && parameters.size() == 2u)
        {
            std::cout << "Test Case 3 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 3 Failed!:\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints4()
    {
        // Concentration and inhibition data for compound input for test case 4 -
        // fitting with
        // more than one data point input
        std::vector<double> testconcentrations = boost::assign::list_of(1.0)(3.0);
        std::vector<double> testinhibitions = boost::assign::list_of(43.35)(70.12);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        // Checks that know parameter results are returned and that exactly 2
        // parameters are returned
        // and prints whether test case has been passed to screen
        if ((fabs(parameters[0] - 1.300) < 1e-3) && (fabs(parameters[1] - 1.020) < 1e-3) && parameters.size() == 2u)
        {
            std::cout << "Test Case 4 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 4 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints5()
    {
        // Concentration and inhibition data for compound input for test case 5 -
        // fitting with
        // more than one data point but all of which have the same concentrations
        std::vector<double> testconcentrations = boost::assign::list_of(10.0)(10.0);
        std::vector<double> testinhibitions = boost::assign::list_of(4.8)(8.1);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        // Checks that know parameter results are returned and that exactly 1
        // parameter is returned
        // and prints whether test case has been passed to screen
        if ((fabs(parameters[0] - 145.039) < 1e-3) && (parameters.size() == 1u))
        {
            std::cout << "Test Case 5 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 5 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints6()
    {
        // Concentration and inhibition data for compound input for test case 5 -
        // fitting with
        // more than one data point but all of which have the same concentrations
        std::vector<double> testconcentrations = boost::assign::list_of(10.0)(10.0)(10.0);
        std::vector<double> testinhibitions = boost::assign::list_of(4.8)(8.1)(6.45);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        // Checks that know parameter results are returned and that exactly 1
        // parameter is returned
        // and prints whether test case has been passed to screen
        if ((fabs(parameters[0] - 145.039) < 1e-3) && (parameters.size() == 1u))

        {
            std::cout << "Test Case 6 Passed!"
                      << "\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 6 Failed!"
                      << "\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints7()
    {
        // Concentration and inhibition data for compound input for test case 5 -
        // fitting with
        // more than one data point but all of which have the same concentrations
        std::vector<double> testconcentrations = boost::assign::list_of(5)(5)(20)(20);
        std::vector<double> testinhibitions = boost::assign::list_of(15)(25)(75)(85);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        // Checks that know parameter results are returned and that exactly 1
        // parameter is returned
        // and prints whether test case has been passed to screen
        if ((fabs(parameters[0] - 10) < 1e-2) && (fabs(parameters[1] - 2) < 1e-2) && (parameters.size() == 2u))

        {
            std::cout << "Test Case 7 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 7 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints8()
    {
        // Concentration and inhibition data for compound input for test case 5 -
        // fitting with
        // more than one data point but all of which have the same concentrations
        std::vector<double> testconcentrations = boost::assign::list_of(5)(5)(5)(20)(20)(20);
        std::vector<double> testinhibitions = boost::assign::list_of(30)(30)(40)(70)(70)(60);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        // Checks that know parameter results are returned and that exactly 1
        // parameter is returned
        // and prints whether test case has been passed to screen
        if ((fabs(parameters[0] - 10) < 1e-2) && (fabs(parameters[1] - 1) < 1e-2) && (parameters.size() == 2u))

        {
            std::cout << "Test Case 8 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 8 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints9()
    {
        std::cout << "Test Case 9 Starting (22 points, with Hill Coeff)... could "
                     "take a minute or two !\n";
        // Wads of data points as per live system - test 1
        std::vector<double> testconcentrations = boost::assign::list_of(
            0.0114311842706904000000)(0.00127013158563227000000)(
            0.308641975308642000000)(0.102880658436214000000)(
            0.000423377195210757000000)(0.0342935528120713000000)(
            0.00381039475689681000000)(0.925925925925926000000)(
            8.33333333333333000000)(2.77777777777778000000)(25.000000)(
            0.0342935528120713000000)(0.00381039475689681000000)(
            2.77777777777778000000)(0.00127013158563227000000)(
            0.0114311842706904000000)(0.308641975308642000000)(
            0.925925925925926000000)(0.000423377195210757000000)(
            8.33333333333333000000)(25.000000)(0.102880658436214000000);
        std::vector<double> testinhibitions = boost::assign::list_of(15)(6)(3)(-1)(-5)(-8)(-10)(-17)(-31)(-32)(-36)(
            11)(10)(8)(3)(2)(-5)(-6)(-6)(-8)(-22)(-44);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        if ((fabs(parameters[0]) > 1000000) && (fabs(parameters[1]) < 5) && (parameters.size() == 2u))
        {
            std::cout << "Test Case 9 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 9 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints10()
    {
        std::cout << "Test Case 10 Starting (almost 100 points, with Hill "
                     "Coeff)... could take a minute or two !\n";
        // Wads of data points as per live system - test 2
        std::vector<double> testconcentrations = boost::assign::list_of(15.000000)(
            0.234375000000)(0.05859375000000)(0.0002288818359375000000)(
            0.003662109375000000)(0.0000143051147460938000000)(0.0146484375000000)(
            0.000057220458984375000000)(0.9375000000)(0.00091552734375000000)(
            3.75000000)(0.0002288818359375000000)(0.003662109375000000)(
            0.000057220458984375000000)(0.00091552734375000000)(0.234375000000)(
            0.05859375000000)(0.0000143051147460938000000)(0.9375000000)(
            0.0146484375000000)(15.000000)(3.75000000)(0.0000143051147460938000000)(
            0.000057220458984375000000)(0.0146484375000000)(0.234375000000)(
            0.05859375000000)(0.00091552734375000000)(0.0002288818359375000000)(
            0.9375000000)(3.75000000)(0.003662109375000000)(15.000000)(
            0.0000143051147460938000000)(0.0002288818359375000000)(
            0.0146484375000000)(0.000057220458984375000000)(0.003662109375000000)(
            0.00091552734375000000)(0.234375000000)(0.05859375000000)(0.9375000000)(
            3.75000000)(15.000000)(0.0146484375000000)(0.003662109375000000)(
            0.234375000000)(0.0000143051147460938000000)(
            0.000057220458984375000000)(0.0002288818359375000000)(0.9375000000)(
            0.00091552734375000000)(0.05859375000000)(15.000000)(3.75000000)(
            0.0146484375000000)(0.0000143051147460938000000)(
            0.00091552734375000000)(0.05859375000000)(0.0002288818359375000000)(
            0.003662109375000000)(0.000057220458984375000000)(3.75000000)(
            0.9375000000)(15.000000)(0.234375000000)(0.00091552734375000000)(
            0.05859375000000)(3.75000000)(0.9375000000)(0.003662109375000000)(
            0.000057220458984375000000)(15.000000)(0.234375000000)(
            0.0146484375000000)(0.0002288818359375000000)(
            0.0000143051147460938000000)(0.0146484375000000)(0.003662109375000000)(
            15.000000)(0.0000143051147460938000000)(0.0002288818359375000000)(
            0.00091552734375000000)(0.000057220458984375000000)(0.234375000000)(
            0.05859375000000)(0.9375000000)(3.75000000)(0.003662109375000000)(
            0.0002288818359375000000)(0.00091552734375000000)(0.234375000000)(
            0.0000143051147460938000000)(15.000000)(0.0146484375000000)(
            0.000057220458984375000000)(0.05859375000000)(0.9375000000)(3.75000000);
        std::vector<double> testinhibitions = boost::assign::list_of(14)(13)(5)(2)(
            -6)(-9)(-9)(-14)(-18)(-25)(-31)(11)(11)(9)(2)(-4)(-5)(-10)(-15)(-17)(
            -19)(-29)(20)(7)(5)(1)(1)(-1)(-3)(-3)(-5)(-12)(-29)(29)(24)(23)(19)(18)(
            9)(5)(3)(-11)(-21)(-24)(35)(28)(19)(17)(15)(13)(10)(10)(-8)(-10)(-14)(
            9)(6)(4)(3)(-3)(-9)(-19)(-19)(-19)(-21)(-38)(10)(8)(4)(3)(0)(-9)(-10)(
            -11)(-22)(-30)(-42)(26)(22)(17)(15)(15)(13)(13)(12)(6)(4)(3)(22)(18)(
            15)(13)(12)(11)(11)(3)(-6)(-8)(-10);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        if ((fabs(parameters[0]) > 1000000) && (fabs(parameters[1]) < 5) && (parameters.size() == 2u))

        {
            std::cout << "Test Case 10 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 10 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints11()
    {
        std::cout << "Test Case 11a Starting (8 points, no Hill Coeff)... could "
                     "take a minute or two !\n";
        std::vector<double> testconcentrations = boost::assign::list_of(
            0.0114311842706904000000)(0.00127013158563227000000)(
            0.308641975308642000000)(0.102880658436214000000)(
            0.000423377195210757000000)(0.0342935528120713000000)(
            0.00381039475689681000000)(0.925925925925926000000);
        std::vector<double> testinhibitions = boost::assign::list_of(15)(6)(3)(-1)(-5)(-8)(-10)(-17);
        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 1u);

        if ((fabs(parameters[0]) > 1000000) && (parameters.size() == 1u))

        {
            std::cout << "Test Case 11a Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 11a Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints11b()
    {
        std::cout << "Test Case 11b Starting (7 points, no Hill Coeff, i.e. as 11a "
                     "data but with one less data point) ... only takes a second "
                     "!\n";
        std::vector<double> testconcentrations = boost::assign::list_of(0.0114311842706904000000)(
            0.00127013158563227000000)(0.308641975308642000000)(
            0.102880658436214000000)(0.000423377195210757000000)(
            0.0342935528120713000000)(0.00381039475689681000000);
        std::vector<double> testinhibitions = boost::assign::list_of(15)(6)(3)(-1)(-5)(-8)(-10);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 1u);

        if ((fabs(parameters[0]) < 16) && (parameters.size() == 1u))

        {
            std::cout << "Test Case 11b Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 11b Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints12()
    {
        std::cout << "Test Case 12 Starting (Ryan's troublesome compounds)\n";
        std::vector<double> testconcentrations = boost::assign::list_of(0.03162)(1)(
            0.1)(100)(0.1)(100)(31.62278)(0.03162)(31.62278)(10)(10)(0.31623)(
            0.03162)(31.62278)(0.1)(0.31623)(0.03162)(3.16228)(100)(3.16228)(
            3.16228)(10)(3.16228)(31.62278)(1)(3.16228)(10)(100)(0.03162)(0.1)(100)(
            100)(10)(10)(0.1)(0.31623)(1)(100)(1)(0.31623)(10)(100)(0.1)(31.62278)(
            1)(3.16228)(0.03162);
        std::vector<double> testinhibitions = boost::assign::list_of(-37.3)(-22.1)(
            -18.7)(-18)(-17.4)(-15.7)(-12.9)(-12.2)(-9.6)(-9)(-7.8)(-7.1)(-6.2)(
            -5.9)(-5.4)(-4.3)(-3.5)(-3.1)(-2.9)(-2.8)(-2.1)(-2)(-1.8)(-1.3)(-1.2)(
            -1)(-0.8)(-0.8)(-0.7)(0)(0)(0.13)(0.47)(0.94)(1.27)(1.77)(2.22)(2.46)(
            2.74)(2.98)(3.05)(3.2)(3.59)(3.9)(6.16)(37.43)(49.03);
        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        if ((fabs(parameters[0]) > 1000000) && (fabs(parameters[1]) > 0.5) && (parameters.size() == 2u))
        {
            std::cout << "Test Case 12 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 12 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints13()
    {
        std::cout << "Test Case 13 Starting (Generated data from a perfect "
                     "dose-response curve)\n";
        std::vector<double> testconcentrations = boost::assign::list_of(0.001)(
            0.003)(0.01)(0.03)(0.1)(0.3)(1)(3)(10)(30)(100)(300);
        std::vector<double> testinhibitions = boost::assign::list_of(
            0.00250808932247271)(0.00968659104658709)(0.0425764008057494)(
            0.164247637072945)(0.718164099678198)(2.71797609381383)(
            10.9404321114155)(32.1792352040277)(67.5975896090632)(88.9597135750537)(
            97.2549035044354)(99.2745211337568);
        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions);

        // Data generated from IC50 of 5.5 and Hill of 1.23
        if ((fabs(parameters[0] - 5.5) < 1e-3) && (fabs(parameters[1] - 1.23) < 1e-4) && (parameters.size() == 2u))
        {
            std::cout << "Test Case 13 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 13 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    // Gef's latest troublesome compound from 7/8/12
    void TestDoseResponseFittingPoints14()
    {
        std::vector<double> testconcentrations = boost::assign::list_of(0.000423377195210757)(0.000423377195210757)(
            0.000423380000000000)(0.00127013000000000)(0.00127013158563227)(
            0.00127013158563227)(0.00381039000000000)(0.00381039475689681)(
            0.00381039475689681)(0.0114311800000000)(0.0114311842706904)(
            0.0114311842706904)(0.0342935500000000)(0.0342935528120713)(
            0.0342935528120713)(0.102880658436214)(0.102880658436214)(
            0.308641975308642)(0.308641975308642)(0.308641980000000)(
            0.925925925925926)(0.925925925925926)(0.925925930000000)(
            2.77777777777778)(2.77777777777778)(2.77777778000000)(
            8.33333333000000)(8.33333333333333)(8.33333333333333)(25)(25)(25);
        //	    std::vector<double> testinhibitions =
        //boost::assign::list_of(-6)(-5)(14.4000000000000)(12.6000000000000)(3)(6)(19.3000000000000)(-10)(10)(18.6000000000000)(2)(15)(11)(-8)(11)(-44)(-1)(-5)(3)(10.5000000000000)(-17)(-6)(14.5000000000000)(-32)(8)(18.1000000000000)(18.2000000000000)(-31)(-8)(-36)(-22)(-0.700000000000000);
        std::vector<double> testinhibitions = boost::assign::list_of(-6.0)(-5)(
            14.4000000000000)(12.6000000000000)(3)(6)(19.3000000000000)(-10)(10)(
            18.6000000000000)(2)(15)(11)(-8)(11)(-44)(-1)(-5)(3)(10.5000000000000)(
            -17)(-6)(14.5000000000000)(-32)(8)(18.1000000000000)(18.2000000000000)(
            -31)(-8)(-36)(-22)(-0.700000000000000);
        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 1u);

        if ((parameters[0] > 1e10) && (parameters.size() == 1u))
        {
            std::cout << "Test Case 14 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 14 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }

        // With IC50 and Hill
        parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 2u);

        if ((parameters[0] > 1e10) && (fabs(parameters[1] - 1.00) < 1e-1) && (parameters.size() == 2u))
        {
            std::cout << "Test Case 14 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 14 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints15()
    {
        std::vector<double> testconcentrations = boost::assign::list_of(
            0.0000143051147460938000000)(0.0002288818359375000000)(
            0.000057220458984375000000)(0.00091552734375000000)(
            0.003662109375000000)(0.0146484375000000)(0.9375000000)(0.234375000000)(
            3.75000000)(0.05859375000000)(15.000000)(0.000057220458984375000000)(
            0.0000143051147460938000000)(0.00091552734375000000)(
            0.0146484375000000)(0.0002288818359375000000)(15.000000)(
            0.003662109375000000)(0.234375000000)(0.9375000000)(0.05859375000000)(
            3.75000000)(15.000000)(3.75000000)(0.234375000000)(0.9375000000)(
            0.0002288818359375000000)(0.00091552734375000000)(
            0.000057220458984375000000)(0.0000143051147460938000000)(
            0.003662109375000000)(0.05859375000000)(0.0146484375000000)(
            0.00091552734375000000)(15.000000)(0.234375000000)(
            0.003662109375000000)(0.0000143051147460938000000)(0.05859375000000)(
            0.000057220458984375000000)(0.9375000000)(3.75000000)(
            0.0002288818359375000000)(0.0146484375000000)(15.000000)(
            0.00091552734375000000)(0.05859375000000)(3.75000000)(
            0.003662109375000000)(0.234375000000)(0.0146484375000000)(
            0.000057220458984375000000)(0.0000143051147460938000000)(0.9375000000)(
            0.0002288818359375000000)(3.75000000)(15.000000)(0.234375000000)(
            0.9375000000)(0.00091552734375000000)(0.0002288818359375000000)(
            0.000057220458984375000000)(0.05859375000000)(0.0146484375000000)(
            0.0000143051147460938000000)(0.003662109375000000);
        std::vector<double> testinhibitions = boost::assign::list_of(21)(21)(20)(
            19)(5)(5)(5)(4)(-1)(-3)(-16)(15)(13)(12)(8)(7)(6)(3)(1)(0)(-1)(-8)(25)(
            15)(12)(11)(10)(10)(9)(2)(1)(1)(-3)(-24)(-24)(-27)(-29)(-29)(-31)(-31)(
            -35)(-37)(-42)(-52)(-2)(-6)(-7)(-8)(-12)(-14)(-19)(-19)(-19)(-20)(-29)(
            4)(3)(0)(-2)(-3)(-7)(-7)(-8)(-8)(-8)(-11);
        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 2u);

        if ((parameters[0] > 1e10) && (parameters.size() == 2u))
        {
            std::cout << "Test Case 15a Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 15a Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }

        // Here we test the rounding of large IC50 values to a sensible cut-off (1e6
        // uM)
        parameters = CallDoseResponseFittingCode(testconcentrations,
                                                 testinhibitions, 2u, true);

        if (fabs(parameters[0] - 1e6) < 1e-6 && fabs(parameters[1] - 1.0) < 1e-6 && (parameters.size() == 2u))
        {
            std::cout << "Test Case 15b Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 15b Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints16()
    {
        std::vector<double> testconcentrations = boost::assign::list_of(
            0.0000143051147460938000000)(0.0002288818359375000000)(
            0.000057220458984375000000)(0.00091552734375000000)(
            0.003662109375000000)(0.0146484375000000)(0.9375000000)(0.234375000000)(
            3.75000000)(0.05859375000000)(15.000000)(0.000057220458984375000000)(
            0.0000143051147460938000000)(0.00091552734375000000)(
            0.0146484375000000)(0.0002288818359375000000)(15.000000)(
            0.003662109375000000)(0.234375000000)(0.9375000000)(0.05859375000000)(
            3.75000000)(15.000000)(3.75000000)(0.234375000000)(0.9375000000)(
            0.0002288818359375000000)(0.00091552734375000000)(
            0.000057220458984375000000)(0.0000143051147460938000000)(
            0.003662109375000000)(0.05859375000000)(0.0146484375000000)(
            0.00091552734375000000)(15.000000)(0.234375000000)(
            0.003662109375000000)(0.0000143051147460938000000)(0.05859375000000)(
            0.000057220458984375000000)(0.9375000000)(3.75000000)(
            0.0002288818359375000000)(0.0146484375000000)(15.000000)(
            0.00091552734375000000)(0.05859375000000)(3.75000000)(
            0.003662109375000000)(0.234375000000)(0.0146484375000000)(
            0.000057220458984375000000)(0.0000143051147460938000000)(0.9375000000)(
            0.0002288818359375000000)(3.75000000)(15.000000)(0.234375000000)(
            0.9375000000)(0.00091552734375000000)(0.0002288818359375000000)(
            0.000057220458984375000000)(0.05859375000000)(0.0146484375000000)(
            0.0000143051147460938000000)(0.003662109375000000);
        std::vector<double> testinhibitions = boost::assign::list_of(21)(21)(20)(
            19)(5)(5)(5)(4)(-1)(-3)(-16)(15)(13)(12)(8)(7)(6)(3)(1)(0)(-1)(-8)(25)(
            15)(12)(11)(10)(10)(9)(2)(1)(1)(-3)(-24)(-24)(-27)(-29)(-29)(-31)(-31)(
            -35)(-37)(-42)(-52)(-2)(-6)(-7)(-8)(-12)(-14)(-19)(-19)(-19)(-20)(-29)(
            4)(3)(0)(-2)(-3)(-7)(-7)(-8)(-8)(-8)(-11);

        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 1u);

        if ((parameters[0] > 1e10) && (parameters.size() == 1u))
        {
            std::cout << "Test Case 16a Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 16a Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }

        // Here we test the rounding of large IC50 values to a sensible cut-off (1e6
        // uM)
        parameters = CallDoseResponseFittingCode(testconcentrations,
                                                 testinhibitions, 1u, true);

        if (fabs(parameters[0] - 1e6) < 1e-6 && (parameters.size() == 1u))
        {
            std::cout << "Test Case 16b Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 16b Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints17a()
    {
        std::vector<double> testconcentrations = boost::assign::list_of(
            15.00000000000000000000)(0.00001430511474609380)(
            3.75000000000000000000)(0.00005722045898437500)(0.23437500000000000000)(
            0.93750000000000000000)(0.00366210937500000000)(0.01464843750000000000)(
            0.00091552734375000000)(0.05859375000000000000)(0.00022888183593750000)(
            0.00091552734375000000)(15.00000000000000000000)(
            0.00366210937500000000)(0.93750000000000000000)(3.75000000000000000000)(
            0.23437500000000000000)(0.00001430511474609380)(0.01464843750000000000)(
            0.00005722045898437500)(0.05859375000000000000)(0.00022888183593750000)(
            0.23437500000000000000)(0.05859375000000000000)(0.00005722045898437500)(
            3.75000000000000000000)(0.93750000000000000000)(0.00022888183593750000)(
            0.00001430511474609380)(0.00366210937500000000)(0.01464843750000000000)(
            15.00000000000000000000)(0.00091552734375000000)(
            0.00366210937500000000)(0.00001430511474609380)(0.00022888183593750000)(
            0.00091552734375000000)(0.93750000000000000000)(0.05859375000000000000)(
            0.01464843750000000000)(0.00005722045898437500)(3.75000000000000000000)(
            15.00000000000000000000)(0.23437500000000000000)(
            0.05859375000000000000)(0.00022888183593750000)(
            15.00000000000000000000)(0.23437500000000000000)(
            0.00366210937500000000)(3.75000000000000000000)(0.93750000000000000000)(
            0.00005722045898437500)(0.01464843750000000000)(0.00091552734375000000)(
            0.00001430511474609380)(0.93750000000000000000)(0.23437500000000000000)(
            0.01464843750000000000)(0.00091552734375000000)(3.75000000000000000000)(
            15.00000000000000000000)(0.00005722045898437500)(
            0.00022888183593750000)(0.00001430511474609380)(0.05859375000000000000)(
            0.00366210937500000000);
        std::vector<double> testinhibitions = boost::assign::list_of(50)(28)(24)(
            21)(19)(15)(15)(14)(13)(13)(4)(28)(22)(22)(20)(20)(19)(19)(16)(10)(4)(
            2)(16)(14)(13)(10)(5)(3)(1)(-9)(-12)(-15)(-20)(34)(23)(21)(21)(20)(19)(
            17)(16)(13)(12)(9)(19)(18)(16)(14)(13)(8)(8)(4)(4)(-3)(-8)(26)(23)(21)(
            19)(17)(15)(-1)(-3)(-5)(-10)(-16);
        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 1u);

        if ((parameters[0] - 54.1465 < 1e-4) && (parameters.size() == 1u))
        {
            std::cout << "Test Case 17a_CaV Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 17a CaV Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints17b()
    {
        std::vector<double> testconcentrations = boost::assign::list_of(
            0.01143118427069040000)(0.00127013158563227000)(0.30864197530864200000)(
            0.10288065843621400000)(0.00042337719521075700)(0.03429355281207130000)(
            0.00381039475689681000)(0.92592592592592600000)(8.33333333333333000000)(
            2.77777777777778000000)(25.00000000000000000000)(
            0.03429355281207130000)(0.00381039475689681000)(2.77777777777778000000)(
            0.00127013158563227000)(0.01143118427069040000)(0.30864197530864200000)(
            0.92592592592592600000)(0.00042337719521075700)(8.33333333333333000000)(
            25.00000000000000000000)(0.10288065843621400000);
        std::vector<double> testinhibitions = boost::assign::list_of(15)(6)(3)(-1)(-5)(-8)(-10)(-17)(-31)(-32)(-36)(
            11)(10)(8)(3)(2)(-5)(-6)(-6)(-8)(-22)(-44);
        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 1u);

        if ((parameters[0] > 1e10) && (parameters.size() == 1u))
        {
            std::cout << "Test Case 17b_KCNQ1 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 17b_KCNQ1 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    // Checking pIC50 rounding
    void TestDoseResponseFittingPoints17b2()
    {
        std::vector<double> testconcentrations = boost::assign::list_of(
            0.01143118427069040000)(0.00127013158563227000)(0.30864197530864200000)(
            0.10288065843621400000)(0.00042337719521075700)(0.03429355281207130000)(
            0.00381039475689681000)(0.92592592592592600000)(8.33333333333333000000)(
            2.77777777777778000000)(25.00000000000000000000)(
            0.03429355281207130000)(0.00381039475689681000)(2.77777777777778000000)(
            0.00127013158563227000)(0.01143118427069040000)(0.30864197530864200000)(
            0.92592592592592600000)(0.00042337719521075700)(8.33333333333333000000)(
            25.00000000000000000000)(0.10288065843621400000);
        std::vector<double> testinhibitions = boost::assign::list_of(15)(6)(3)(-1)(-5)(-8)(-10)(-17)(-31)(-32)(-36)(
            11)(10)(8)(3)(2)(-5)(-6)(-6)(-8)(-22)(-44);
        std::vector<double> parameters = CallDoseResponseFittingCode(
            testconcentrations, testinhibitions, 1u, true);

        if ((parameters[0] = 1e6) && (parameters.size() == 1u))
        {
            std::cout << "Test Case 17b2_KCNQ1 Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 17b2_KCNQ1 Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingPoints17c()
    {
        std::vector<double> testconcentrations = boost::assign::list_of(
            11.11111111111110000000)(0.13717421124828500000)(
            0.04572473708276180000)(0.00508052634252909000)(0.00169350878084303000)(
            3.70370370370370000000)(1.23456790123457000000)(0.41152263374485600000)(
            33.33333333333330000000)(0.01524157902758730000)(
            100.00000000000000000000)(33.33333333333330000000)(
            1.23456790123457000000)(0.04572473708276180000)(
            11.11111111111110000000)(100.00000000000000000000)(
            0.13717421124828500000)(0.41152263374485600000)(3.70370370370370000000)(
            0.00508052634252909000)(0.01524157902758730000)(0.00169350878084303000);
        std::vector<double> testinhibitions = boost::assign::list_of(24)(22)(18)(
            14)(14)(11)(11)(10)(10)(9)(3)(22)(16)(14)(13)(13)(12)(10)(7)(7)(4)(3);
        std::vector<double> parameters = CallDoseResponseFittingCode(testconcentrations, testinhibitions, 1u);

        if ((parameters[0] - 546.424 < 1e-4) && (parameters.size() == 1u))
        {
            std::cout << "Test Case 17c_NaV Passed!\n";
            PrintParameterValues(parameters);
        }
        else
        {
            std::cout << "Test Case 17c_NaV Failed!\n";
            PrintParameterValues(parameters);
            TS_ASSERT(false);
        }
    }

    void TestDoseResponseFittingWithLimits()
    {
        std::vector<double> doses = boost::assign::list_of(1)(10);
        std::vector<double> responses = boost::assign::list_of(10)(20);

        unsigned num_params_to_fit = 2u;
        bool round_massive_ic50s = true;

        RunHillFunctionMinimization testcompound(
            doses, responses, num_params_to_fit, round_massive_ic50s);
        std::vector<double> params = testcompound.Run();

        TS_ASSERT_EQUALS(params.size(), 2u);
        TS_ASSERT_DELTA(params[0], 512.2845, 1e-4);
        TS_ASSERT_DELTA(params[1], 0.3521, 1e-4);

        testcompound.SetHillLimits(0.6, 5);

        params = testcompound.Run();

        TS_ASSERT_EQUALS(params.size(), 2u);
        TS_ASSERT_DELTA(params[0], 88.3678, 1e-4);
        TS_ASSERT_DELTA(params[1], 0.6000, 1e-4);
    }
};
#endif /*TESTDOSERESPONSEFITTING_HPP_*/
