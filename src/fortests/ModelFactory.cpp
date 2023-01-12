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

#include "ModelFactory.hpp"

ModelFactory::TModelMapping ModelFactory::mpModelRegistry;

ModelFactory::TModelMapping ModelFactory::getModelRegistry()
{
    if (ModelFactory::mpModelRegistry == nullptr)
    {
        ModelFactory::mpModelRegistry = std::make_unique<std::map<std::pair<std::string, std::string>, ModelFactory::TCreateMethod>>();
    }
    return ModelFactory::mpModelRegistry;
}

bool ModelFactory::Exists(const std::string& rName, const std::string& rType)
{
    std::pair<std::string, std::string> name_type = std::make_pair(rName, rType);
    return ModelFactory::getModelRegistry()->find(name_type) != ModelFactory::getModelRegistry()->end();
}

void* ModelFactory::Create(const std::string& rName, 
                           const std::string& rType, 
                           boost::shared_ptr<AbstractIvpOdeSolver> pSolver, 
                           boost::shared_ptr<AbstractStimulusFunction> pStimulus)
{
    std::pair<std::string, std::string> name_type = std::make_pair(rName, rType);
    auto it = ModelFactory::getModelRegistry()->find(name_type);
    if (it == ModelFactory::getModelRegistry()->end())
    {
        EXCEPTION("Model type combination does not exist cannot create: " + rName + ", " + rType);
    }
    return it->second(pSolver, pStimulus); // call the createFunc
}

bool ModelFactory::Register(const std::string& rName, 
                            const std::string& rType, 
                            ModelFactory::TCreateMethod funcCreate)
{
    std::pair<std::string, std::string> nameType = std::make_pair(rName, rType);
    if (ModelFactory::getModelRegistry()->count(nameType) != 0)
    {
        EXCEPTION("Duplicate model: " + rName + " registration with the ModelFactory for type: " + rType + ". If you are using your own version of this model please rename the cellml file.");
    }
    ModelFactory::getModelRegistry()->insert(std::pair<std::pair<std::string, std::string>, ModelFactory::TCreateMethod>(nameType, funcCreate));
    return true;
}
