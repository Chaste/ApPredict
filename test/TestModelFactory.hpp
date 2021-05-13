/*

Copyright (c) 2005-2021, University of Oxford.
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

#ifdef CHASTE_CVODE

#ifndef _TESTMODELFACTORY_HPP_
#define _TESTMODELFACTORY_HPP_

#include "CommandLineArgumentsMocker.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "SimpleStimulus.hpp"

#include "ModelFactory.hpp"
#include "SetupModel.hpp"

#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "mahajan_shiferaw_2008Cvode.hpp"
#include "hund_rudy_2004Cvode.hpp"
#include "grandi_pasqualini_bers_2010_ssCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "paci_hyttinen_aaltosetala_severi_ventricularVersionCvode.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"

#include "ohara_rudy_cipa_v1_2017.hpp"
#include "ohara_rudy_cipa_v1_2017BackwardEulerOpt.hpp"

class TestModelFactory : public CxxTest::TestSuite
{
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    std::vector<std::unique_ptr<AbstractCvodeCell>> referenceModels;

    void checkModels(std::vector<boost::shared_ptr<AbstractCvodeCell>> models) {
	for(unsigned int i=0; i < models.size(); i++){
	    TS_ASSERT(models[i] != nullptr);
	    TS_ASSERT(models[i]->GetNumberOfStateVariables() == referenceModels[i]->GetNumberOfStateVariables());
	    TS_ASSERT(models[i]->GetVoltage() == referenceModels[i]->GetVoltage());
	    TS_ASSERT(models[i]->GetNumberOfParameters() == referenceModels[i]->GetNumberOfParameters());
	}
    }


public:
    void setUp() {
	p_stimulus.reset(new SimpleStimulus(-25.5, 2.0, 50.0));
	p_solver.reset(new EulerIvpOdeSolver);

	referenceModels.push_back(std::make_unique<Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode>(p_solver, p_stimulus));
	referenceModels.push_back(std::make_unique<Cellten_tusscher_model_2006_epiFromCellMLCvode>(p_solver, p_stimulus));
	referenceModels.push_back(std::make_unique<Cellmahajan_shiferaw_2008FromCellMLCvode>(p_solver, p_stimulus));
	referenceModels.push_back(std::make_unique<Cellhund_rudy_2004FromCellMLCvode>(p_solver, p_stimulus));
	referenceModels.push_back(std::make_unique<Cellgrandi_pasqualini_bers_2010_ssFromCellMLCvode>(p_solver, p_stimulus));
	referenceModels.push_back(std::make_unique<Cellohara_rudy_2011_endoFromCellMLCvode>(p_solver, p_stimulus));
	referenceModels.push_back(std::make_unique<Cellpaci_hyttinen_aaltosetala_severi_ventricularVersionFromCellMLCvode>(p_solver, p_stimulus));
	referenceModels.push_back(std::make_unique<Cellohara_rudy_cipa_v1_2017FromCellMLCvode>(p_solver, p_stimulus));
     }

    void TestModelFactoryName()
    {
	boost::shared_ptr<AbstractCvodeCell> model1((AbstractCvodeCell*)ModelFactory::Create("shannon_wang_puglisi_weber_bers_2004" , "AnalyticCvode", p_solver, p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> model2((AbstractCvodeCell*)ModelFactory::Create("ten_tusscher_model_2006_epi" , "AnalyticCvode", p_solver, p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> model3((AbstractCvodeCell*)ModelFactory::Create("mahajan_shiferaw_2008" , "AnalyticCvode", p_solver, p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> model4((AbstractCvodeCell*)ModelFactory::Create("hund_rudy_2004" , "AnalyticCvode", p_solver, p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> model5((AbstractCvodeCell*)ModelFactory::Create("grandi_pasqualini_bers_2010_ss" , "AnalyticCvode", p_solver, p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> model6((AbstractCvodeCell*)ModelFactory::Create("ohara_rudy_2011_endo" , "AnalyticCvode", p_solver, p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> model7((AbstractCvodeCell*)ModelFactory::Create("paci_hyttinen_aaltosetala_severi_ventricularVersion" , "AnalyticCvode", p_solver, p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> model8((AbstractCvodeCell*)ModelFactory::Create("ohara_rudy_cipa_v1_2017" , "AnalyticCvode", p_solver, p_stimulus));

        checkModels({model1, model2, model3, model4, model5, model6, model7, model8});
    }

    void TestWrongModelOrType()
    {
	boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
	boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

        std::unique_ptr<AbstractCvodeCell> model1((AbstractCvodeCell*)ModelFactory::Create("wrong_model_name" , "AnalyticCvode", p_solver, p_stimulus));
        std::unique_ptr<AbstractCvodeCell> model2((AbstractCvodeCell*)ModelFactory::Create("shannon_wang_puglisi_weber_bers_2004" , "Wrong model type", p_solver, p_stimulus));

        TS_ASSERT(model1 == nullptr);
        TS_ASSERT(model2 == nullptr);
    }

    void TestModelFactorySetupModel()
    {
        SetupModel setup1(1.0, 1u);
        SetupModel setup2(1.0, 2u);
        SetupModel setup3(1.0, 3u);
        SetupModel setup4(1.0, 4u);
        SetupModel setup5(1.0, 5u);
        SetupModel setup6(1.0, 6u);
        SetupModel setup7(1.0, 7u);
        SetupModel setup8(1.0, 8u);

        checkModels({setup1.GetModel(), setup2.GetModel(), setup3.GetModel(), setup4.GetModel(), setup5.GetModel(), setup6.GetModel(), setup7.GetModel(), setup8.GetModel()});
    }

    void TestSetupModelCommandLineNumber()
    {
        CommandLineArgumentsMocker wrapper("--model 1");
        SetupModel setup1(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper2("--model 2");
        SetupModel setup2(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper3("--model 3");
        SetupModel setup3(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper4("--model 4");
        SetupModel setup4(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper5("--model 5");
        SetupModel setup5(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper6("--model 6");
        SetupModel setup6(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper7("--model 7");
        SetupModel setup7(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper8("--model 8");
        SetupModel setup8(1.0, UNSIGNED_UNSET);

        checkModels({setup1.GetModel(), setup2.GetModel(), setup3.GetModel(), setup4.GetModel(), setup5.GetModel(), setup6.GetModel(), setup7.GetModel(), setup8.GetModel()});
    }

    void TestSetupModelCommandLineName()
    {
        CommandLineArgumentsMocker wrapper("--model shannon_wang_puglisi_weber_bers_2004");
        SetupModel setup1(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper2("--model ten_tusscher_model_2006_epi");
        SetupModel setup2(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper3("--model mahajan_shiferaw_2008");
        SetupModel setup3(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper4("--model hund_rudy_2004");
        SetupModel setup4(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper5("--model grandi_pasqualini_bers_2010_ss");
        SetupModel setup5(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper6("--model ohara_rudy_2011_endo");
        SetupModel setup6(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper7("--model paci_hyttinen_aaltosetala_severi_ventricularVersion");
        SetupModel setup7(1.0, UNSIGNED_UNSET);
        CommandLineArgumentsMocker wrapper8("--model ohara_rudy_cipa_v1_2017");
        SetupModel setup8(1.0, UNSIGNED_UNSET);

        checkModels({setup1.GetModel(), setup2.GetModel(), setup3.GetModel(), setup4.GetModel(), setup5.GetModel(), setup6.GetModel(), setup7.GetModel(), setup8.GetModel()});
    }

    void TestSetupModelCommandLineWrongIndex()
    {
        CommandLineArgumentsMocker wrapper2("--model 99999");
        TS_ASSERT_THROWS_THIS(SetupModel setup1(1.0, UNSIGNED_UNSET),
                              "No model matches this index: 99999");

        CommandLineArgumentsMocker wrapper("--model unknownmodel");
        TS_ASSERT_THROWS_THIS(SetupModel setup1(1.0, UNSIGNED_UNSET),
                              "No model matches this index: unknownmodel");
    }

    void TestReRegisterExistingModel()
    {
        ModelFactory::Register("TestModel", "TestType", nullptr);
        TS_ASSERT_THROWS_THIS(ModelFactory::Register("TestModel", "TestType", nullptr),
                              "Duplicate model: TestModel registration with the ModelFactory for type: TestType. If you are using your own version of this model please rename the cellml file.");
    }


    void TestOtherModelTypes()
    {
        // Check we are actually generating the correct type of the models
        // Set up Normal, Cvode and BE models, compute 1 ms and compare
        Cellohara_rudy_cipa_v1_2017FromCellML normal(p_solver, p_stimulus);
        Cellohara_rudy_cipa_v1_2017FromCellMLCvode cvode(p_solver, p_stimulus);
        Cellohara_rudy_cipa_v1_2017FromCellMLBackwardEulerOpt be(p_solver, p_stimulus);

        // Create Normal, Cvode and BE models via factory
        boost::shared_ptr<AbstractCardiacCell> factoryNormal((AbstractCardiacCell*)ModelFactory::Create("ohara_rudy_cipa_v1_2017" , "Normal", p_solver, p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> factoryCvode((AbstractCvodeCell*)ModelFactory::Create("ohara_rudy_cipa_v1_2017" , "AnalyticCvode", p_solver, p_stimulus));
        boost::shared_ptr<AbstractBackwardEulerCardiacCell<21>> factoryBe((AbstractBackwardEulerCardiacCell<21>*)ModelFactory::Create("ohara_rudy_cipa_v1_2017" , "BackwardEulerOpt", p_solver, p_stimulus));

        // Compute 1 ms
        normal.Compute(0.0, 1.0);
        cvode.Compute(0.0, 1.0);
        be.Compute(0.0, 1.0);
        factoryNormal->Compute(0.0, 1.0);
        factoryCvode->Compute(0.0, 1.0);
        factoryBe->Compute(0.0, 1.0);

        // The different model types and their Factory created version should have the same voltage
        TS_ASSERT(normal.GetVoltage() == factoryNormal->GetVoltage());
        TS_ASSERT(cvode.GetVoltage() == factoryCvode->GetVoltage());
        TS_ASSERT(be.GetVoltage() == factoryBe->GetVoltage());

        // The voltage for the different types should be very close but not identical
        TS_ASSERT_DELTA(factoryNormal->GetVoltage(), factoryCvode->GetVoltage(), 1e-4);
        TS_ASSERT_DELTA(factoryNormal->GetVoltage(), factoryBe->GetVoltage(), 1e-4);

        TS_ASSERT(factoryCvode->GetVoltage() != factoryNormal->GetVoltage());
        TS_ASSERT(factoryBe->GetVoltage() != factoryNormal->GetVoltage());
        TS_ASSERT(factoryCvode->GetVoltage() != factoryBe->GetVoltage());
    }

};

#endif //_TESTMODELFACTORY_HPP_
#endif //_CHASTE_CVODE
