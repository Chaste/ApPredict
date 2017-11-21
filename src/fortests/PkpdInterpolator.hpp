/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef PKPDINTERPOLATOR_HPP_
#define PKPDINTERPOLATOR_HPP_

#include <boost/shared_ptr.hpp>
#include "PkpdDataStructure.hpp"

/**
 * A class to run some simulations for PKPD, and interpolate values back at concentrations
 * in the PKPD file. 
 */
class PkpdInterpolator
{
  private:
    friend class TestPkpdInterpolations; // To test interpolations.

    /** A data reader class to hold information read from the PKPD file */
    boost::shared_ptr<PkpdDataStructure> mpPkpdReader;
    
    /** 
      * Perform linear interpolation to get an estimate of y_star at x_star 
      * @param x_star The independent variable at which to get an interpolated value
      * @param rX  The vector of independent variables.
      * @param rY  The vector of dependent variables to interpolate between.
      */
    double DoLinearInterpolation(double x_star,const std::vector<double>& rX, const std::vector<double>& rY) const;

  public:
    /** Constructor */
    PkpdInterpolator();
    
    /** The main method to perform all the analysis for this class, 
     * Runs ApPredictMethods and then does interpolation on the results 
     */
    void Run();
};

#endif // PKPDINTERPOLATOR_HPP
