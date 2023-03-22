# Copyright (c) 2005-2023, University of Oxford.
# All rights reserved.
# 
# University of Oxford means the Chancellor, Masters and Scholars of the
# University of Oxford, having an administrative office at Wellington
# Square, Oxford OX1 2JD, UK.
# 
# This file is part of Chaste.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * Neither the name of the University of Oxford nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# Example CMake build configuration file for a Chaste user project.
# See also https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/CmakeBuildGuide for full documentation.

# You can set which Chaste components (or other projects) your project depends on by editing the
# find_package() call for Chaste.  E.g.
#   find_package(Chaste COMPONENTS cell_based)
# for a project just using the cell_based component (and its dependencies), or
#   find_package(Chaste COMPONENTS heart lung)
# for a project combining heart & lung simulations.
#
# Note that the order in which components are specified does not matter.

include(FetchContent)
#message (STATUS "Downloading APPredict cellml files to ${CMAKE_CURRENT_SOURCE_DIR}/src/cellml")
FetchContent_Declare(CellML_repo_appredict
                     GIT_REPOSITORY https://github.com/Chaste/cellml.git
                     GIT_TAG c26e1a0939e8c619e2ba10fe41d1af495503f6a1
                     SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/cellml)
FetchContent_MakeAvailable(CellML_repo_appredict)

# remove .git folder to prevent accidentally checking in deletion of cellml files
file(REMOVE_RECURSE ${CMAKE_CURRENT_SOURCE_DIR}/src/cellml/.git)

# remove cellml files that are not required
file(GLOB cellml_files ${CMAKE_CURRENT_SOURCE_DIR}/src/cellml/cellml/*.cellml)
foreach(cellml_file ${cellml_files})
    if (NOT ${cellml_file} MATCHES davies_isap_2012.cellml
        AND NOT ${cellml_file} MATCHES decker_2009.cellml
        AND NOT ${cellml_file} MATCHES faber_rudy_2000.cellml
        AND NOT ${cellml_file} MATCHES grandi_pasqualini_bers_2010_ss.cellml
        AND NOT ${cellml_file} MATCHES hund_rudy_2004.cellml
        AND NOT ${cellml_file} MATCHES livshitz_rudy_2007.cellml
        AND NOT ${cellml_file} MATCHES mahajan_shiferaw_2008.cellml
        AND NOT ${cellml_file} MATCHES ohara_rudy_2011_endo.cellml
        AND NOT ${cellml_file} MATCHES ohara_rudy_cipa_v1_2017.cellml
        AND NOT ${cellml_file} MATCHES paci_hyttinen_aaltosetala_severi_ventricularVersion.cellml
        AND NOT ${cellml_file} MATCHES priebe_beuckelmann_1998.cellml
        AND NOT ${cellml_file} MATCHES shannon_wang_puglisi_weber_bers_2004.cellml
        AND NOT ${cellml_file} MATCHES ten_tusscher_model_2006_epi.cellml)
        file(REMOVE ${cellml_file})
    endif()
endforeach()

# Here we just depend on core components (nothing application-specific).
find_package(Chaste COMPONENTS heart)

# Alternatively, to specify a Chaste installation directory use a line like that below.
# This is needed if your project is not contained in the projects folder within a Chaste source tree.
#find_package(Chaste COMPONENTS heart crypt PATHS /path/to/chaste-install NO_DEFAULT_PATH)

# Here we add extra arguments to force PyCML to use this extra argument (make Get and Set methods for
# all metadata annotated variables). 
set(Chaste_CODEGEN_EXTRA_ARGS "--use-model-factory")

# ApPredict uses pthreads, so we specify that threads are required and add the necessary compile flag
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
if(CMAKE_USE_PTHREADS_INIT)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

# Change the project name in the line below to match the folder this file is in,
# i.e. the name of your project.
chaste_do_project(ApPredict)
