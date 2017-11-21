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

void PkpdInterpolator::Run()
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

    std::vector<double> concs;
    std::vector<double> apd90s;

    // Restrict the scope of ApPredictMethods to avoid any funny file conflicts etc.
    // and ensure it is tidy before doing PKPD.
    {
        ApPredictMethods ap_predict;
        ap_predict.SetMaxConcentrationForPkpd(max_conc);
        ap_predict.Run();

        concs = ap_predict.GetConcentrations();
        apd90s = ap_predict.GetApd90s();
    }

    std::cout << "Conc\tApd90" << std::endl;
    for (unsigned i = 0; i < concs.size(); i++)
    {
        std::cout << concs[i] << "\t" << apd90s[i] << std::endl;
    }

    // Open up the existing output folder
    OutputFileHandler handler("ApPredict_output", false); // false = don't wipe

    // Copy the PKPD data file into the output folder for posterity
    FileFinder pkpd_file(CommandLineArguments::Instance()->GetStringCorrespondingToOption("--pkpd-file"),
                         RelativeTo::AbsoluteOrCwd);
    handler.CopyFileTo(pkpd_file);

    // Open a results file - in a try catch as it is conceivable someone could have named their PK file this!
    out_stream p_output_file;
    try
    {
        p_output_file = handler.OpenOutputFile("pkpd_results.txt");
    }
    catch (Exception& e)
    {
        EXCEPTION("Could not open a new output file called pkpd_results.txt. Error was " << e.GetMessage());
    }
    *p_output_file << "Time";
    for (unsigned i = 0; i < mpPkpdReader->GetNumberOfPatients(); i++)
    {
        *p_output_file << "\tConc_for_patient_" << i << "(uM)";
    }
    *p_output_file << std::endl;

    std::vector<double> times = mpPkpdReader->GetTimes();
    for (unsigned i = 0; i < times.size(); i++)
    {
        *p_output_file << times[i];
        const std::vector<double>& r_concs_at_this_time = mpPkpdReader->GetConcentrationsAtTimeIndex(i);
        for (unsigned p=0; p<r_concs_at_this_time.size(); p++)
        {
            double interpolated_apd90 = DoLinearInterpolation(r_concs_at_this_time[p],concs,apd90s);
            *p_output_file << "\t" << interpolated_apd90;
        }
        *p_output_file << std::endl;
    }
    p_output_file->close();
}

/* Perform linear interpolation to get an estimate of y_star at x_star */
double PkpdInterpolator::DoLinearInterpolation(double x_star, const std::vector<double>& rX, const std::vector<double>& rY) const
{
    if (x_star <= rX[0])
    {
        return rY[0];
    }
    if (x_star >= rX.back())
    {
        return rY.back();
    }
    
    auto lower = std::lower_bound(rX.cbegin(), rX.cend(), x_star);
    unsigned lower_idx = lower - rX.cbegin();
    
    // I'll let the compiler tidy all this up!
    double lower_x = rX[lower_idx-1u];
    double upper_x = rX[lower_idx];
    double lower_y = rY[lower_idx-1u];
    double upper_y = rY[lower_idx];

    return lower_y + ((x_star-lower_x)/(upper_x-lower_x))*(upper_y-lower_y);
}


