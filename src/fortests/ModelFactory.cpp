#include "ModelFactory.hpp"
#include <iostream>

std::unique_ptr<std::map<std::pair<std::string, std::string>, ModelFactory::TCreateMethod>> modelRegistry = std::make_unique<std::map<std::pair<std::string, std::string>, ModelFactory::TCreateMethod>>();

void* ModelFactory::Create(const std::string& name, const std::string& type, boost::shared_ptr<AbstractIvpOdeSolver> pSolver, boost::shared_ptr<AbstractStimulusFunction> pStimulus)
{
    std::pair<std::string, std::string> name_type = std::make_pair(name, type);
    auto it = modelRegistry->find(name_type);
    if (it != modelRegistry->end())
    {
        return it->second(pSolver, pStimulus); // call the createFunc
    }

    return nullptr;
}

bool ModelFactory::Register(const std::string& name, const std::string& type, ModelFactory::TCreateMethod funcCreate)
{
    std::pair<std::string, std::string> nameType = std::make_pair(name, type);
    if (modelRegistry->count(nameType) == 0)
    {
        modelRegistry->insert(std::pair<std::pair<std::string, std::string>, ModelFactory::TCreateMethod>(nameType, funcCreate));
        return true;
    }
    return false;
}