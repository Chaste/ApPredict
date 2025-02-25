/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef MODELFACTORY_HPP_
#define MODELFACTORY_HPP_

#include <memory>
#include <string>
#include <map>
#include "AbstractIvpOdeSolver.hpp"
#include "RegularStimulus.hpp"


/**
 * Class to provide a way to create models indirectly using the name of the cellml 
 * file and the model type. This eliminated the need to hardcode which models are available.
 * Following the Factory design pattern, this class provides methods for generated 
 * cells to register their existance and methods for creating new instances of cells.
 */
class ModelFactory
{
public:
    /** Type definition for methods to creat a new model */
    using TCreateMethod = void*(*)(boost::shared_ptr<AbstractIvpOdeSolver> p_solver, boost::shared_ptr<AbstractStimulusFunction> p_stimulus);

    /** Type definition for Mapping of model name & type to create method */
    using TModelMapping = std::shared_ptr<std::map<std::pair<std::string, std::string>, TCreateMethod>>;

    /** Method to check whether a given model, type combination has been registered 
     *  and can be created */
    static bool Exists(const std::string& rName, const std::string& rType);

    /** Method to create an given type instance of a model. 
     * The result needs to be type cast to the desired cell type 
     * (see TestModelFactory.hpp for examples) */
    static void* Create(const std::string& rName, 
                        const std::string& rType, 
                        boost::shared_ptr<AbstractIvpOdeSolver> pSolver, 
                        boost::shared_ptr<AbstractStimulusFunction> pStimulus);

    /** Registration method for generated cells to register themselves with the model factory */
    static bool Register(const std::string& rName, 
                         const std::string& rType, 
                         ModelFactory::TCreateMethod funcCreate);

private:
    // N.B. these need to be down here as the types are declared above.

    /** Mapping of model name & type to create method */
    static ModelFactory::TModelMapping mpModelRegistry;

    /** Private getter, making sure the map is initialised before using it */
    static TModelMapping  getModelRegistry();
};

#endif // MODELFACTORY_HPP_
