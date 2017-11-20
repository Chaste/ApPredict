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

#include "PkpdInterpolator.hpp"

#include "ApPredictMethods.hpp"
#include "CommandLineArguments.hpp"
#include "Exception.hpp"

PkpdInterpolator::PkpdInterpolator()
{
    if (!CommandLineArguments::Instance()->OptionExists("--pkpd-file"))
    {
        EXCEPTION("PkpdInterpolator class needs a PKPD file to be specified with --pkpd-file <file_path> argument.");
    }

    FileFinder pkpd_file(CommandLineArguments::Instance()->GetStringCorrespondingToOption("--pkpd-file"),
                         RelativeTo::AbsoluteOrCwd);

    if (!pkpd_file.IsFile())
    {
        EXCEPTION("The File '" << pkpd_file.GetAbsolutePath() << "' does not exist. Please give a relative or absolute path.");
    }

    mpPkpdReader = boost::shared_ptr<PkpdDataStructure>(new PkpdDataStructure(pkpd_file));
}

void PkpdInterpolator::RunApPredict()
{
    // Calculate maximum concentration to use... add 10% to the maximum we saw.
    double max_conc = 1.1 * mpPkpdReader->GetMaximumConcentration();

    std::cout << "Max conc to use = " << max_conc << " uM" << std::endl;

    if (CommandLineArguments::Instance()->OptionExists("--plasma-conc-high"))
    {
        EXCEPTION("The argument --plasma-conc-high will be ignored. Using PKPD file to set concentrations. Please remove it to avoid confustion!");
    }

    if (CommandLineArguments::Instance()->OptionExists("--plasma-concs"))
    {
        EXCEPTION("The argument --plasma-concs will be ignored. Using PKPD file to set concentrations. Please remove it to avoid confustion!");
    }

    ApPredictMethods ap_predict;
    ap_predict.SetMaxConcentrationForPkpd(max_conc);
    ap_predict.Run();

    std::vector<double> concs = ap_predict.GetConcentrations();
    std::vector<double> apd90s = ap_predict.GetApd90s();

    std::cout << "Conc\tApd90" << std::endl;
    for (unsigned i = 0; i < concs.size(); i++)
    {
        std::cout << concs[i] << "\t" << apd90s[i] << std::endl;
    }
}
