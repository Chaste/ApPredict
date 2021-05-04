#ifndef MODELFACTORY_HPP_
#define MODELFACTORY_HPP_

#include <memory>
#include <string>
#include <map>
#include "AbstractIvpOdeSolver.hpp"
#include "RegularStimulus.hpp"


class ModelFactory{
public:
    using TCreateMethod = void*(*)(boost::shared_ptr<AbstractIvpOdeSolver> p_solver, boost::shared_ptr<AbstractStimulusFunction> p_stimulus);
    static void* Create(const std::string& name, const std::string& type, boost::shared_ptr<AbstractIvpOdeSolver> p_solver, boost::shared_ptr<AbstractStimulusFunction> p_stimulus);
    static bool Register(const std::string name, const std::string& type, TCreateMethod funcCreate);
};
#endif // MODELFACTORY_HPP_