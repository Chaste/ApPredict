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

#include"NelderMeadMinimizer.hpp"
#include"HillFunction.hpp"
#include <cmath>
#include <vector>
#include<iostream>

//Constructor for class with default variables
NelderMeadMinimizer::NelderMeadMinimizer(std::vector<double>& rParameters, HillFunction* p_hill_function)
: 	mrParameters(rParameters),
  	mpFunctionToMinimise(p_hill_function),
    mNumParameters(rParameters.size()),
    mMaxNumIterations(100000000u),
    mReflectionCoefficient(1.0),
    mTolerance(1e-8),
    mContractionCoefficient(0.5),
    mExpansionCoefficient(2.0),
    mShrinkCoefficient(0.5),
    mDisplayIterations(false),
    mNumFunctionEvaluations(0u)
{
	// Defining default variables above that are used in the Nelder-Mead algorithm
}

//Method in which minimization algorithm using the Nelder-Mead simplex method is implemented
void NelderMeadMinimizer::Minimize()
{
	mNumFunctionEvaluations = 0;

	//Declares variables and the vertices matrix used in the algorithm
	double **vertices; //Matrix storing vertices of simplex
	unsigned smallest_value_vertex=0; //index of vertex with smallest error value (just initialised to please compiler)
	unsigned next_smallest_value_vertex; // index of vertex with next smallest error value
	unsigned largest_value_vertex; // index of vertex with largest error value

	double function_reflection_value; //Value of error related to function at reflection point
	double function_expansion_value; //Value of error related to function at expansion point
	double function_contraction_value; //Value of error related to function at contraction point

	double centroidsum; //Variable for sum of individual terms for sum for calculating centroid
	//double functionsum; //Variables for testing for convergence - for old stopping criteria
	//double functionaverage;
	//double averagesquarediff;

	//Dynamically allocating vertex matrix of simplex
	vertices = new double* [mNumParameters+1]; // The simplex is a 'tetrahedron-ish' in R^n, so n+1 vertices
	for (unsigned i=0; i<mNumParameters+1; i++)
	{
		vertices[i] = new double [mNumParameters]; // Each vertex is defined by n parameter values.
	}

	//Allocating standard vector for the function value vector - that is error associated with fitting hill function at each vertex point
	std::vector<double> function_value(mNumParameters+1);

	//Allocating standard vector for the reflection_vertex vector - updated vertex points following reflection
	std::vector<double> reflection_vertex(mNumParameters);

	//Allocating standard vector for the expansion_vertex vector - updated vertex points following expansion
	std::vector<double> expansion_vertex(mNumParameters);

	//Allocating standard vector for the contraction_vertex vector - updated vertex points following contraction
	std::vector<double> contraction_vertex(mNumParameters);

	//Allocating standard vector for the centroid_vertex vector
	std::vector<double> centroid_vertex(mNumParameters);

	//Vector of initial guess at solution
	std::vector<double> x_0(mNumParameters);
	x_0[0] = mrParameters[0];
	if (mNumParameters ==2)
	{
		x_0[1] = mrParameters[1];
	}

	// Pfeffer method - as used in fminsearch - handles
	//cases such as x_0[0] = 1 and x_0[1] = 10^5 much better
	//than some other alternative initial simplex formulations

	double delta_u = 0.05;
	double delta_z = 0.00025;

	//Create the initial simplex, assuming one of the vertices is (0,0)

	//(0,0) vertex
	for (unsigned i=0; i<mNumParameters; i++)
	{
		vertices[0][i] = x_0[i];
	}
	// Initialisation of other vertices
	for (unsigned i=1; i<=mNumParameters;i++)
	{
		for (unsigned j=0; j<mNumParameters;j++)
		{
			if (i-1==j&&x_0[j]!=0)
			{
				vertices[i][j]=(1+delta_u)*x_0[j];
			}
			if (i-1==j && x_0[j]==0)
			{
				vertices[i][j]=delta_z;
			}

			if (i-1!=j)
			{
				vertices[i][j] = x_0[j];
			}
		}
	}

	// Beginning of loop over iterations
	for (unsigned iteration=0; iteration<mMaxNumIterations; iteration++)
	{
		if (iteration>0 && mDisplayIterations)
		{
			std::cout << "Iteration = " << iteration << ", FuncEvals = " << mNumFunctionEvaluations << ", min f(x) =  " << function_value[smallest_value_vertex] << "\n" <<std::flush;
		}
		// Now finding the initial function values - that is the minimum value error of the hill
		// function fit and expected value at each vertex point.
		// Appropriate points of the vertices matrix are put into a vector (argument) to pass to evaluation of
		// associated error when fitting to hill function

		for (unsigned j=0; j<=mNumParameters; j++)
		{
			std::vector<double> argument;
			for (unsigned i=0; i<mNumParameters; i++)
			{
				argument.push_back(vertices[j][i]);
			}
			function_value[j] = mpFunctionToMinimise->Evaluate(argument);
			mNumFunctionEvaluations++;
		}
		// Ordering the vertices so that they are ordered in order of ascending value of the error given
		// after fitting to the hill function

		//Finding vertex associated with the largest error

		largest_value_vertex = 0;
		for (unsigned j=0; j<=mNumParameters; j++)
		{
			if (function_value[j]>function_value[largest_value_vertex])
			{
				largest_value_vertex = j;
			}
		}

		// Finding vertex associated with smallest error
		smallest_value_vertex = 0u;
		for (unsigned j=0; j<=mNumParameters; j++)
		{
			if (function_value[j] < function_value[smallest_value_vertex])
			{
				smallest_value_vertex = j;
			}
		}

		// Ordering remaining points according to error calculated at each vertex
		next_smallest_value_vertex = smallest_value_vertex;

		for (unsigned j=0; j<=mNumParameters; j++)
		{
			if (   (function_value[j] > function_value[next_smallest_value_vertex])
					&& (function_value[j] < function_value[largest_value_vertex])       )
			{
				next_smallest_value_vertex = j;
			}
		}

		// Calculating the centroid of the n vertices which give the least error

		for(unsigned j=0; j<=mNumParameters-1; j++)
		{
			centroidsum=0.0;

			for(unsigned m=0; m<=mNumParameters;m++)
			{
				if (m!=largest_value_vertex)
				{
					centroidsum += vertices[m][j];

				}
			}
			centroid_vertex[j] = centroidsum/mNumParameters;
		}

		// Computing the reflection point by reflecting largest_value_vertex to the new vertex reflection_vertex

		for (unsigned j=0; j<=mNumParameters-1;j++)
		{
			reflection_vertex[j] = (1 + mReflectionCoefficient)*centroid_vertex[j]
			                                                                    - mReflectionCoefficient*vertices[largest_value_vertex][j];
		}

		// Sending to evaluate method to calculate the error associated value of the hill function with the
		// parameters at the reflection vertex and updating function reflection value
		function_reflection_value = mpFunctionToMinimise->Evaluate(reflection_vertex);
		mNumFunctionEvaluations++;

		// Calculating expansion vertex and evaluating function value at the point
		if (function_reflection_value < function_value[smallest_value_vertex])
		{
			for (unsigned j=0; j<=mNumParameters-1;j++)
			{
				expansion_vertex[j] = (1+mReflectionCoefficient*mExpansionCoefficient)*centroid_vertex[j]
				                                                                                       -(mReflectionCoefficient*mExpansionCoefficient)*vertices[largest_value_vertex][j];
			}
			function_expansion_value = mpFunctionToMinimise->Evaluate(expansion_vertex);
			mNumFunctionEvaluations++;

			//Terminate iteration and accept expansion point if criteria satisfied
			if (function_expansion_value<function_reflection_value)
			{
				for (unsigned j=0;j<mNumParameters;j++)
				{
					vertices[largest_value_vertex][j] = expansion_vertex[j];

				}
				function_value[largest_value_vertex] = function_expansion_value;
				//	std::cout<<"Expansion"<<"\n";
				goto terminate;

			}
			// Otherwise accept reflection point and terminate iteration
			else
			{
				for (unsigned j=0;j<mNumParameters;j++)
				{
					vertices[largest_value_vertex][j] = reflection_vertex[j];
				}
				function_value[largest_value_vertex] = function_reflection_value;
				//	std::cout<<"Reflection"<<"\n";
				goto terminate;
			}
		}
		//if function_value[smallest_value_vertex]<=function_reflection_value and appropriate
		//criteria satisfied then accept the reflected point and terminate the iteration
		else
		{
			if (function_reflection_value<function_value[next_smallest_value_vertex])
			{
				for (unsigned j=0;j<mNumParameters;j++)
				{
					vertices[largest_value_vertex][j] = reflection_vertex[j];
				}
				function_value[largest_value_vertex] = function_reflection_value;
				//	std::cout<<"Reflection"<<"\n";
				goto terminate;
			}
			else // if function_reflection_value>=function_value[next_smallest_value_vertex]
				// then perform a contraction
			{
				// Perform outside contraction if this criteria is satisfied
				if (function_reflection_value<function_value[largest_value_vertex])
				{
					for (unsigned j=0;j<mNumParameters;j++)
					{
						contraction_vertex[j] = (1+mReflectionCoefficient*mContractionCoefficient)*centroid_vertex[j]
						                                                                                           - (mContractionCoefficient*mReflectionCoefficient)*vertices[largest_value_vertex][j];
					}

					function_contraction_value = mpFunctionToMinimise->Evaluate(contraction_vertex);
					mNumFunctionEvaluations++;

					// If this criteria satisfied accept outside contraction point and terminate iteration
					if (function_contraction_value<=function_reflection_value)
					{
						for (unsigned j=0;j<mNumParameters;j++)
						{
							vertices[largest_value_vertex][j] = contraction_vertex[j];
						}
						function_value[largest_value_vertex] = function_contraction_value;
						//	std::cout<<"Outside Contraction"<<"\n";
						goto terminate;
					}

					else //Perform a shrink step
					{

						for (unsigned j=0;j<=mNumParameters-1;j++)
						{
							vertices[next_smallest_value_vertex][j] =  vertices[smallest_value_vertex][j]
							                                                                           + mShrinkCoefficient*(vertices[next_smallest_value_vertex][j]- vertices[smallest_value_vertex][j]);
							vertices[largest_value_vertex][j] =  vertices[smallest_value_vertex][j] + mShrinkCoefficient*(vertices[largest_value_vertex][j]- vertices[smallest_value_vertex][j]);


						}

						//	std::cout<<"Shrink"<<"\n";

					}
				}
				else //Perform an inside contraction
				{
					for (unsigned j=0;j<mNumParameters;j++)
					{
						contraction_vertex[j] = (1-mContractionCoefficient)*centroid_vertex[j]
						                                                                    + mContractionCoefficient*vertices[largest_value_vertex][j];
					}

					function_contraction_value = mpFunctionToMinimise->Evaluate(contraction_vertex);
					mNumFunctionEvaluations++;

					//If criteria satisfied then accept inside contraction point and terminate iteration
					if(function_value[largest_value_vertex]>function_contraction_value)
					{
						for(unsigned j=0; j<mNumParameters; j++)
						{
							vertices[largest_value_vertex][j] = contraction_vertex[j];
						}
						function_value[largest_value_vertex] = function_contraction_value;
						//	std::cout<<"Inside Contraction"<<"\n";

						goto terminate;
					}

					else //Perform a shrink step
					{
						for (unsigned j=0;j<=mNumParameters-1;j++)
						{
							vertices[next_smallest_value_vertex][j] =  vertices[smallest_value_vertex][j] + mShrinkCoefficient*(vertices[next_smallest_value_vertex][j]- vertices[smallest_value_vertex][j]);
							vertices[largest_value_vertex][j] =  vertices[smallest_value_vertex][j] + mShrinkCoefficient*(vertices[largest_value_vertex][j]- vertices[smallest_value_vertex][j]);

						}

						//	std::cout<<"Shrink"<<"\n";

					}
				}


			}


		}

		// Relevant vertices matrix points are put into the argument vector to be passed to
		// Evaluate for evaluation using the hill function and vertices are reordered appropriately
		// following shrink step

		{
			std::vector<double> argument;
			for (unsigned i=0; i<mNumParameters; i++)
			{
				argument.push_back(vertices[next_smallest_value_vertex][i]);
			}
			function_value[next_smallest_value_vertex] = mpFunctionToMinimise->Evaluate(argument);
			mNumFunctionEvaluations++;
		}

		{
			std::vector<double> argument;
			for (unsigned i=0; i<mNumParameters; i++)
			{
				argument.push_back(vertices[largest_value_vertex][i]);
			}
			function_value[largest_value_vertex] = mpFunctionToMinimise->Evaluate(argument);
			mNumFunctionEvaluations++;
		}
		/*		// Testing for convergence of minimal error between hill function value and expected value associated with optimised parameter solution
		// Stopping criterion as proposed by Nelder & Mead
		functionsum = 0.0;

		for (unsigned j=0;j<=mNumParameters;j++)
		{
			functionsum +=function_value[j];
		}

		functionaverage = functionsum/(mNumParameters+1);

		averagesquarediff = 0.0;
		for (unsigned j=0; j<=mNumParameters;j++)
		{
			averagesquarediff += pow((function_value[j]-functionaverage),2.0)/(mNumParameters);

		}
		averagesquarediff = sqrt(averagesquarediff);

		// Terminates once the required tolerance is satisfied
		if (averagesquarediff<mTolerance)
		{
			//std::cout<<"Number of Iterations = "<< iteration <<"\n";
			break;
		}*/

		terminate:
		// Stopping criteria as used in fminsearch

		// Establish maximum difference in vertex values from smallest value vertex - to
		// be used in stopping criteria
		double v = -10;
		double max_vertices_difference;
		max_vertices_difference = 10000;
		for (unsigned i = 0; i<mNumParameters; i++)
		{
			if ((fabs(vertices[largest_value_vertex][i]-vertices[smallest_value_vertex][i]))>v)
			{
				max_vertices_difference = fabs(vertices[largest_value_vertex][i]-vertices[smallest_value_vertex][i]);
			}
			else
			{
				max_vertices_difference = v;
			}

			v = max_vertices_difference;

		}

		for (unsigned i = 0; i<mNumParameters; i++)
		{
			if ((fabs(vertices[next_smallest_value_vertex][i]-vertices[smallest_value_vertex][i]))>v)
			{
				max_vertices_difference = fabs(vertices[next_smallest_value_vertex][i]-vertices[smallest_value_vertex][i]);
			}
			else
			{
				max_vertices_difference = v;
			}

			v = max_vertices_difference;
		}

		//Stopping criteria
		if(((fabs(function_value[largest_value_vertex]-function_value[smallest_value_vertex])<=mTolerance)&&
				(fabs(function_value[next_smallest_value_vertex]-function_value[smallest_value_vertex])<=mTolerance))&&
				(max_vertices_difference <= mTolerance))

		{
			break;
		}
	}

	if (mDisplayIterations)	std::cout << "Simplex minimisation complete\n";

	// end of minimisation loop

	// Parameters updated once minimisation complete
	for (unsigned i=0; i<mNumParameters;i++)
	{
		mrParameters[i] = vertices[smallest_value_vertex][i];
		if (mDisplayIterations) std::cout<< "Param[" << i << "] = " << mrParameters[i]<<"\n";
	}
	// Deleting dynamically allocated memory
	for (unsigned i=0;i<mNumParameters+1;i++)
	{
		delete[] vertices[i];
	}
	delete[] vertices;
}

