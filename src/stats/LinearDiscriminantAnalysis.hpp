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

#ifndef LINEARDISCRIMINANTANALYSIS_HPP_
#define LINEARDISCRIMINANTANALYSIS_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Exception.hpp"

using namespace boost::numeric::ublas;

class LinearDiscriminantAnalysis
{
private:
    /** Needed for serialization. */
    friend class TestLinearDiscriminantAnalysis;

    /**
     * A standard vector of matrices of training data for each group
     *
     * rows of matrices correspond to training points,
     * columns correspond to discriminant variables.
     */
    const std::vector<matrix<double> > mTraining;

    /** The dimension of the classification */
    unsigned mDimension;

    std::vector<vector<double> > mMeanTrainingPoints;

    std::vector<matrix<double> > mCovarianceMatrices;

    matrix<double> mPooledCovarianceMatrix;

    std::vector<vector<double> > mInvPooledDotMean;

    /**
     * Matrix inversion routine.
     * Uses lu_factorize and lu_substitute in uBLAS to invert a matrix.
     *
     * This is based upon Numerical Recipes in C...
     *
     * @param rInput  The matrix to be inverted
     * @param rInverse  The matrix to be populated with the inverse.
     */
    void InvertMatrix(const matrix<double>& rInput, matrix<double>& rInverse);

    /**
     * Get the Covariance Matrix of a given input matrix
     *
     * @param rInput  The matrix to use
     * @param rInverse  The matrix to be populated with the covariant (must be the correct size already).
     */
    void CalculateCovariance(const matrix<double>& rInput, matrix<double>& rCov);

    /**
     * @param rTraining  The training data, each entry in the vector is a different group.
     * @param rCovMats  Covariance matrices for each group (empty - filled by this method)
     * @param rPooledCov  The pooled covariance matrix (empty - filled by this method).
     */
    void CalculatePooledCovariance(const std::vector<matrix<double> >& rTraining,
                          std::vector<matrix<double> >& rCovMats,
                          matrix<double>& rPooledCov);

    /**
     * Calculate the mean point in each of the training groups.
     *
     * @param rTraining  The training data (mTraining is put in here by internal methods).
     * @param rMeans  Gets filled up by this method
     */
    void CalculateMeanPoints(const std::vector<matrix<double> >& rTraining, std::vector<vector<double> >& rMeans);

public:
    /**
     * Normal constructor
     *
     * @param rTraining  a vector of training data. Each entry corresponds to one of the training groups.
     * @param testing  Whether we are conducting testing on the methods and do not wish to run constructor properly.
     */
    LinearDiscriminantAnalysis(const std::vector<matrix<double> >& rTraining, bool testing=false);

    /**
     * @return  The mean point in each of the training groups.
     */
    std::vector<vector<double> > GetMeanTrainingPoints();

    /**
     * @return  The pooled covariance matrix of the training data.
     */
    matrix<double> GetPooledCovarianceMatrix();

    /**
     * @return  The covariance matrices of each training group.
     */
    std::vector<matrix<double> > GetCovarianceMatrices();

    /**
     * Perform Linear Discriminant Analysis.
     *
     * @param rPoint  The point to classify.
     * @return The index of the training group to which rPoint has been assigned.
     */
    unsigned ClassifyThisPoint(const vector<double>& rPoint);

};

#endif /* LINEARDISCRIMINANTANALYSIS_HPP_ */
