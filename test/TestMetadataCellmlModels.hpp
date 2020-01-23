/*

Copyright (c) 2005-2020, University of Oxford.
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

#ifndef _TESTMETADATACELLMLMODELS_HPP_
#define _TESTMETADATACELLMLMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include "Exception.hpp"

#include <boost/shared_ptr.hpp>
#include "AbstractCvodeCell.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ZeroStimulus.hpp"

#include "Shannon2004.hpp"
#include "davies_isap_2012.hpp"
#include "davies_isap_2012Cvode.hpp"
#include "davies_isap_2012CvodeOpt.hpp"
#include "davies_isap_2012Opt.hpp"
#include "hund_rudy_2004.hpp"
#include "hund_rudy_2004Cvode.hpp"
#include "hund_rudy_2004CvodeOpt.hpp"
#include "hund_rudy_2004Opt.hpp"
#include "livshitz_rudy_2007.hpp"
#include "livshitz_rudy_2007Cvode.hpp"
#include "livshitz_rudy_2007CvodeOpt.hpp"
#include "livshitz_rudy_2007Opt.hpp"
#include "mahajan_shiferaw_2008.hpp"
#include "mahajan_shiferaw_2008Cvode.hpp"
#include "mahajan_shiferaw_2008CvodeOpt.hpp"
#include "mahajan_shiferaw_2008Opt.hpp"
#include "priebe_beuckelmann_1998.hpp"
#include "priebe_beuckelmann_1998Cvode.hpp"
#include "priebe_beuckelmann_1998CvodeOpt.hpp"
#include "priebe_beuckelmann_1998Opt.hpp"
//#include "Shannon2004Opt.hpp" // Opt files not made to save time as Shannon is
//in trunk.
#include "Shannon2004Cvode.hpp"
//#include "Shannon2004CvodeOpt.hpp" // Opt files not made to save time as
//Shannon is in trunk.
#include "decker_2009.hpp"
#include "decker_2009Cvode.hpp"
#include "decker_2009CvodeOpt.hpp"
#include "decker_2009Opt.hpp"
#include "grandi_pasqualini_bers_2010_ss.hpp"
#include "grandi_pasqualini_bers_2010_ssCvode.hpp"
#include "grandi_pasqualini_bers_2010_ssCvodeOpt.hpp"
#include "grandi_pasqualini_bers_2010_ssOpt.hpp"
#include "ohara_rudy_2011_endo.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "ohara_rudy_2011_endoCvodeOpt.hpp"
#include "ohara_rudy_2011_endoOpt.hpp"
#include "ohara_rudy_cipa_v1_2017.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017CvodeOpt.hpp"
#include "ohara_rudy_cipa_v1_2017Opt.hpp"
#include "paci_hyttinen_aaltosetala_severi_ventricularVersion.hpp"
#include "paci_hyttinen_aaltosetala_severi_ventricularVersionCvode.hpp"
#include "paci_hyttinen_aaltosetala_severi_ventricularVersionCvodeOpt.hpp"
#include "paci_hyttinen_aaltosetala_severi_ventricularVersionOpt.hpp"
#include "ten_tusscher_model_2006_epi.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_epiCvodeOpt.hpp"
#include "ten_tusscher_model_2006_epiOpt.hpp"

class TestMetadataCellmlModels : public CxxTest::TestSuite
{
public:
    /**
	* This test just ensures that all cell models can be compiled and given
	* metadata values.
	*
	* (i.e. that pyCML and the metadata tags generated valid Chaste cell models.)
	*/
    void TestMetadataHasBeenCorrectlyTranslatedToChaste(void)
    {
#ifdef CHASTE_CVODE
        // Regular stimulus and an ode solver object
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);

        // We use the Get/SetParameter() methods now to intervene in the right hand
        // side of the ODE.
        double a = 0.9;
        double b = 0.8;
        double c = 0.7;
        double d = 0.6;
        double e = 0.5;
        double f = 0.4;

        for (unsigned model_index = 1u; model_index <= 12u; model_index++)
        {
            AbstractCardiacCell* p_chaste_cell;
            AbstractCardiacCell* p_chaste_cell_opt;
            AbstractCvodeCell* p_cvode_cell;
            AbstractCvodeCell* p_cvode_cell_opt;

            std::cout << "Model " << model_index << ": ";
            switch (model_index)
            {
                case 1u:
                {
                    p_chaste_cell = new CellShannon2004FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new CellShannon2004FromCellML(
                        p_solver,
                        p_stimulus); // Opt files not made as Shannon is in trunk
                    p_cvode_cell = new CellShannon2004FromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new CellShannon2004FromCellMLCvode(
                        p_solver,
                        p_stimulus); // Opt files not made as Shannon is in trunk
                    break;
                }
                case 2u:
                {
                    p_chaste_cell = new Cellten_tusscher_model_2006_epiFromCellML(
                        p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellten_tusscher_model_2006_epiFromCellMLOpt(
                        p_solver, p_stimulus);
                    p_cvode_cell = new Cellten_tusscher_model_2006_epiFromCellMLCvode(
                        p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellten_tusscher_model_2006_epiFromCellMLCvodeOpt(p_solver,
                                                                                             p_stimulus);
                    break;
                }
                case 3u:
                {
                    p_chaste_cell = new Cellpriebe_beuckelmann_1998FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellpriebe_beuckelmann_1998FromCellMLOpt(
                        p_solver, p_stimulus);
                    p_cvode_cell = new Cellpriebe_beuckelmann_1998FromCellMLCvode(
                        p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellpriebe_beuckelmann_1998FromCellMLCvodeOpt(
                        p_solver, p_stimulus);
                    break;
                }
                case 4u:
                {
                    p_chaste_cell = new Cellmahajan_shiferaw_2008FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellmahajan_shiferaw_2008FromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Cellmahajan_shiferaw_2008FromCellMLCvode(
                        p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellmahajan_shiferaw_2008FromCellMLCvodeOpt(
                        p_solver, p_stimulus);
                    break;
                }
                case 5u:
                {
                    p_chaste_cell = new Celllivshitz_rudy_2007FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Celllivshitz_rudy_2007FromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Celllivshitz_rudy_2007FromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Celllivshitz_rudy_2007FromCellMLCvodeOpt(
                        p_solver, p_stimulus);
                    break;
                }
                case 6u:
                {
                    p_chaste_cell = new Cellhund_rudy_2004FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellhund_rudy_2004FromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Cellhund_rudy_2004FromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellhund_rudy_2004FromCellMLCvodeOpt(p_solver, p_stimulus);
                    break;
                }
                case 7u:
                {
                    p_chaste_cell = new Cellgrandi_pasqualini_bers_2010_ssFromCellML(
                        p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellgrandi_pasqualini_bers_2010_ssFromCellMLOpt(p_solver,
                                                                                            p_stimulus);
                    p_cvode_cell = new Cellgrandi_pasqualini_bers_2010_ssFromCellMLCvode(
                        p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellgrandi_pasqualini_bers_2010_ssFromCellMLCvodeOpt(
                        p_solver, p_stimulus);
                    break;
                }
                case 8u:
                {
                    p_chaste_cell = new Celldecker_2009FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Celldecker_2009FromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Celldecker_2009FromCellMLCvodeOpt(p_solver, p_stimulus);
                    break;
                }
                case 9u:
                {
                    p_chaste_cell = new Cellohara_rudy_2011_endoFromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellohara_rudy_2011_endoFromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellohara_rudy_2011_endoFromCellMLCvodeOpt(
                        p_solver, p_stimulus);
                    break;
                }
                case 10u:
                {
                    p_chaste_cell = new Cellohara_rudy_cipa_v1_2017FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellohara_rudy_cipa_v1_2017FromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Cellohara_rudy_cipa_v1_2017FromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellohara_rudy_cipa_v1_2017FromCellMLCvodeOpt(
                        p_solver, p_stimulus);
                    break;
                }
                case 11u:
                {
                    p_chaste_cell = new Celldavies_isap_2012FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Celldavies_isap_2012FromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Celldavies_isap_2012FromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Celldavies_isap_2012FromCellMLCvodeOpt(p_solver, p_stimulus);
                    break;
                }
                case 12u:
                {
                    p_chaste_cell = new Cellpaci_hyttinen_aaltosetala_severi_ventricularVersionFromCellML(
                        p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellpaci_hyttinen_aaltosetala_severi_ventricularVersionFromCellMLOpt(
                        p_solver, p_stimulus);
                    p_cvode_cell = new Cellpaci_hyttinen_aaltosetala_severi_ventricularVersionFromCellMLCvode(
                        p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellpaci_hyttinen_aaltosetala_severi_ventricularVersionFromCellMLCvodeOpt(
                        p_solver, p_stimulus);
                    break;
                }
                default:
                {
                    EXCEPTION("Unrecognised cell model");
                }
            }
            std::cout << p_chaste_cell->GetSystemName() << std::endl
                      << std::flush;
            TS_ASSERT_EQUALS(p_chaste_cell->HasParameter(
                                 "membrane_fast_sodium_current_conductance"),
                             true);
            TS_ASSERT_EQUALS(p_chaste_cell->HasParameter(
                                 "membrane_L_type_calcium_current_conductance"),
                             true);
            TS_ASSERT_EQUALS(
                p_chaste_cell->HasParameter(
                    "membrane_rapid_delayed_rectifier_potassium_current_conductance"),
                true);
            TS_ASSERT_EQUALS(
                p_chaste_cell->HasParameter(
                    "membrane_slow_delayed_rectifier_potassium_current_conductance"),
                true);

            const double default_a = p_chaste_cell->GetParameter(
                "membrane_fast_sodium_current_conductance");
            const double default_b = p_chaste_cell->GetParameter(
                "membrane_L_type_calcium_current_conductance");
            const double default_c = p_chaste_cell->GetParameter(
                "membrane_rapid_delayed_rectifier_potassium_current_conductance");
            const double default_d = p_chaste_cell->GetParameter(
                "membrane_slow_delayed_rectifier_potassium_current_conductance");

            double default_e = 0;
            if (model_index == 7u)
            {
                TS_ASSERT_EQUALS(
                    p_chaste_cell->HasParameter(
                        "membrane_fast_transient_outward_current_conductance"),
                    true);
                default_e = p_chaste_cell->GetParameter(
                    "membrane_fast_transient_outward_current_conductance");
            }
            else if (model_index != 5u)
            {
                TS_ASSERT_EQUALS(p_chaste_cell->HasParameter(
                                     "membrane_transient_outward_current_conductance"),
                                 true);
                default_e = p_chaste_cell->GetParameter(
                    "membrane_transient_outward_current_conductance");
            }

            double default_f = 0;
            if (model_index == 9u || model_index == 10u) // O'Hara and variant
            {
                default_f = p_chaste_cell->GetParameter("membrane_persistent_sodium_current_conductance");
            }

            {
                p_chaste_cell->SetParameter("membrane_fast_sodium_current_conductance",
                                            default_a * a);
                p_chaste_cell->SetParameter(
                    "membrane_L_type_calcium_current_conductance", default_b * b);
                p_chaste_cell->SetParameter(
                    "membrane_rapid_delayed_rectifier_potassium_current_conductance",
                    default_c * c);
                p_chaste_cell->SetParameter(
                    "membrane_slow_delayed_rectifier_potassium_current_conductance",
                    default_d * d);
                if (model_index == 7u)
                {
                    p_chaste_cell->SetParameter(
                        "membrane_fast_transient_outward_current_conductance",
                        default_e * e);
                }
                else if (model_index != 5u)
                {
                    p_chaste_cell->SetParameter(
                        "membrane_transient_outward_current_conductance", default_e * e);
                }

                if (model_index == 9u || model_index == 10u) // O'Hara and variant
                {
                    p_chaste_cell->SetParameter("membrane_persistent_sodium_current_conductance", f * default_f);
                }

                TS_ASSERT_EQUALS(p_chaste_cell->HasCellMLDefaultStimulus(), true);
                p_chaste_cell->UseCellMLDefaultStimulus();
            }

            {
                p_chaste_cell_opt->SetParameter(
                    "membrane_fast_sodium_current_conductance", default_a * a);
                p_chaste_cell_opt->SetParameter(
                    "membrane_L_type_calcium_current_conductance", default_b * b);
                p_chaste_cell_opt->SetParameter(
                    "membrane_rapid_delayed_rectifier_potassium_current_conductance",
                    default_c * c);
                p_chaste_cell_opt->SetParameter(
                    "membrane_slow_delayed_rectifier_potassium_current_conductance",
                    default_d * d);
                if (model_index == 7u)
                {
                    p_chaste_cell_opt->SetParameter(
                        "membrane_fast_transient_outward_current_conductance",
                        default_e * e);
                }
                else if (model_index != 5u)
                {
                    p_chaste_cell_opt->SetParameter(
                        "membrane_transient_outward_current_conductance", default_e * e);
                }

                if (model_index == 9u || model_index == 10u) // O'Hara and variant
                {
                    p_chaste_cell_opt->SetParameter("membrane_persistent_sodium_current_conductance", f * default_f);
                }

                TS_ASSERT_EQUALS(p_chaste_cell_opt->HasCellMLDefaultStimulus(), true);
                p_chaste_cell_opt->UseCellMLDefaultStimulus();
            }

            {
                p_cvode_cell->SetParameter("membrane_fast_sodium_current_conductance",
                                           default_a * a);
                p_cvode_cell->SetParameter(
                    "membrane_L_type_calcium_current_conductance", default_b * b);
                p_cvode_cell->SetParameter(
                    "membrane_rapid_delayed_rectifier_potassium_current_conductance",
                    default_c * c);
                p_cvode_cell->SetParameter(
                    "membrane_slow_delayed_rectifier_potassium_current_conductance",
                    default_d * d);
                if (model_index == 7u)
                {
                    p_cvode_cell->SetParameter(
                        "membrane_fast_transient_outward_current_conductance",
                        default_e * e);
                }
                else if (model_index != 5u)
                {
                    p_cvode_cell->SetParameter(
                        "membrane_transient_outward_current_conductance", default_e * e);
                }

                if (model_index == 9u || model_index == 10u) // O'Hara and variant
                {
                    p_cvode_cell->SetParameter("membrane_persistent_sodium_current_conductance", f * default_f);
                }

                TS_ASSERT_EQUALS(p_cvode_cell->HasCellMLDefaultStimulus(), true);
                p_cvode_cell->UseCellMLDefaultStimulus();
            }

            {
                p_cvode_cell_opt->SetParameter(
                    "membrane_fast_sodium_current_conductance", default_a * a);
                p_cvode_cell_opt->SetParameter(
                    "membrane_L_type_calcium_current_conductance", default_b * b);
                p_cvode_cell_opt->SetParameter(
                    "membrane_rapid_delayed_rectifier_potassium_current_conductance",
                    default_c * c);
                p_cvode_cell_opt->SetParameter(
                    "membrane_slow_delayed_rectifier_potassium_current_conductance",
                    default_d * d);
                if (model_index == 7u)
                {
                    p_cvode_cell_opt->SetParameter(
                        "membrane_fast_transient_outward_current_conductance",
                        default_e * e);
                }
                else if (model_index != 5u)
                {
                    p_cvode_cell_opt->SetParameter(
                        "membrane_transient_outward_current_conductance", default_e * e);
                }
                if (model_index == 9u || model_index == 10u) // O'Hara and variant
                {
                    p_cvode_cell_opt->SetParameter("membrane_persistent_sodium_current_conductance", f * default_f);
                }

                TS_ASSERT_EQUALS(p_cvode_cell_opt->HasCellMLDefaultStimulus(), true);
                p_cvode_cell_opt->UseCellMLDefaultStimulus();
            }

            delete p_chaste_cell;
            delete p_chaste_cell_opt;
            delete p_cvode_cell;
            delete p_cvode_cell_opt;
        }
#else
        std::cout << "You need CVODE installed and Chaste hostconfig set up to use "
                     "it for this project.\n";
#endif
    }
};

#endif //_TESTMETADATACELLMLMODELS_HPP_
