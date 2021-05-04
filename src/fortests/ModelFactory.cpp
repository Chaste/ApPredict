#include "ModelFactory.hpp"
#include <iostream>

std::unique_ptr<std::map<std::pair<std::string, std::string>, ModelFactory::TCreateMethod>> s_methods = std::make_unique<std::map<std::pair<std::string, std::string>, ModelFactory::TCreateMethod>>();

void* ModelFactory::Create(const std::string& name, const std::string& type, boost::shared_ptr<AbstractIvpOdeSolver> p_solver, boost::shared_ptr<AbstractStimulusFunction> p_stimulus)
{
    std::pair<std::string, std::string> nameType = std::make_pair(name, type);
    auto it = s_methods->find(nameType);
    if (it != s_methods->end())
    {
        return it->second(p_solver, p_stimulus); // call the createFunc
    }

    return nullptr;
}

bool ModelFactory::Register(const std::string name, const std::string& type, ModelFactory::TCreateMethod funcCreate)
{
    std::pair<std::string, std::string> nameType = std::make_pair(name, type);
    if (s_methods->count(nameType) == 0)
    {
        std::cout << name << " " << funcCreate << std::endl;
        s_methods->insert(std::pair<std::pair<std::string, std::string>, ModelFactory::TCreateMethod>(nameType, funcCreate));
        return true;
    }
    return false;
}