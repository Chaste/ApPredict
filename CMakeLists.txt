# Copyright (c) 2005-2025, University of Oxford.
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
# See also https://chaste.github.io/docs/dev-guides/cmake-build-guide for full documentation.

# You can set which Chaste components (or other projects) your project depends on by editing the
# find_package() call for Chaste.  E.g.
#   find_package(Chaste COMPONENTS cell_based)
# for a project just using the cell_based component (and its dependencies), or
#   find_package(Chaste COMPONENTS heart lung)
# for a project combining heart & lung simulations.
#
# Note that the order in which components are specified does not matter.

set(ApPredict_CELLML_TAG fb76f77ef973ff4e18c2b36b841ff3be5e59b5c8)
set(ApPredict_CELLML_ARCHIVE_URL "https://github.com/Chaste/cellml/archive/${ApPredict_CELLML_TAG}.tar.gz")
set(ApPredict_CELLML_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/cellml")
set(ApPredict_CELLML_FILES
    davies_isap_2012.cellml
    decker_2009.cellml
    faber_rudy_2000.cellml
    grandi_pasqualini_bers_2010_ss.cellml
    hund_rudy_2004.cellml
    livshitz_rudy_2007.cellml
    mahajan_shiferaw_2008.cellml
    ohara_rudy_2011_endo.cellml
    ohara_rudy_cipa_v1_2017.cellml
    paci_hyttinen_aaltosetala_severi_ventricularVersion.cellml
    priebe_beuckelmann_1998.cellml
    shannon_wang_puglisi_weber_bers_2004.cellml
    ten_tusscher_model_2006_epi.cellml
    noble_model_1962.cellml
)

# Retry cellml downloads a few times with a short backoff before giving up.
set(ApPredict_CELLML_DOWNLOAD_ATTEMPTS 3)

# Work out which cellml files need downloading.
set(cellml_missing)
foreach(cellml_file ${ApPredict_CELLML_FILES})
    if (NOT EXISTS "${ApPredict_CELLML_SRC}/${cellml_file}")
        list (APPEND cellml_missing ${cellml_file})
    endif ()
endforeach()

if (cellml_missing)
    # Fetch the cellml repository as a tarball rather than fetching individual files
    # (the total download size is smaller, and there is less chance of rate limiting).
    set(cellml_tarball "${CMAKE_CURRENT_BINARY_DIR}/cellml-${ApPredict_CELLML_TAG}.tar.gz")

    message (STATUS "Downloading ApPredict cellml archive ${ApPredict_CELLML_ARCHIVE_URL}")

    set(dl_attempt 0)
    set(dl_status_code 1)
    while (NOT (dl_status_code EQUAL 0) AND (dl_attempt LESS ApPredict_CELLML_DOWNLOAD_ATTEMPTS))
        math(EXPR dl_attempt "${dl_attempt} + 1")
        if (dl_attempt GREATER 1)
            execute_process(COMMAND "${CMAKE_COMMAND}" -E sleep ${dl_attempt})
        endif ()
        file(DOWNLOAD
            "${ApPredict_CELLML_ARCHIVE_URL}"
            "${cellml_tarball}"
            TIMEOUT 60
            TLS_VERIFY ON
            STATUS dl_status
        )
        list (GET dl_status 0 dl_status_code)
        list (GET dl_status 1 dl_status_msg)
        if (NOT dl_status_code EQUAL 0)
            file(REMOVE "${cellml_tarball}") # remove stub
            message (STATUS "Attempt ${dl_attempt} of ${ApPredict_CELLML_DOWNLOAD_ATTEMPTS} to download the cellml archive failed: ${dl_status_msg}")
        endif ()
    endwhile ()

    # Failure to get the cellml files is fatal as the build cannot complete without them.
    if (NOT dl_status_code EQUAL 0)
        message (FATAL_ERROR "Could not download ${ApPredict_CELLML_ARCHIVE_URL} after ${dl_attempt} attempts: ${dl_status_msg}")
    endif ()

    # Unpack the downloaded cellml archive.
    set(cellml_tmp "${CMAKE_CURRENT_BINARY_DIR}/cellml-unpack")
    file(REMOVE_RECURSE "${cellml_tmp}")
    file(MAKE_DIRECTORY "${cellml_tmp}")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E tar xzf "${cellml_tarball}"
        WORKING_DIRECTORY "${cellml_tmp}"
        RESULT_VARIABLE untar_result
    )
    if (NOT untar_result EQUAL 0)
        message (FATAL_ERROR "Could not unpack the cellml archive ${cellml_tarball}")
    endif ()

    # Copy the required cellml files into the source tree.
    message (STATUS "Copying ApPredict cellml files to ${ApPredict_CELLML_SRC}")
    foreach(cellml_file ${cellml_missing})
        # Anchored on the cellml directory so a similarly-named model cannot match.
        file(GLOB cellml_found "${cellml_tmp}/*/cellml/${cellml_file}")
        if (NOT cellml_found)
            message (FATAL_ERROR "${cellml_file} is not present in the cellml repository at revision ${ApPredict_CELLML_TAG}")
        endif ()
        list (GET cellml_found 0 cellml_found_file)
        file(COPY "${cellml_found_file}" DESTINATION "${ApPredict_CELLML_SRC}")
    endforeach()

    # Clean up.
    file(REMOVE_RECURSE "${cellml_tmp}")
    file(REMOVE "${cellml_tarball}")
endif ()

# Here we just depend on core components (nothing application-specific).
find_package(Chaste COMPONENTS heart)

# Alternatively, to specify a Chaste installation directory use a line like that below.
# This is needed if your project is not contained in the projects folder within a Chaste source tree.
#find_package(Chaste COMPONENTS heart crypt PATHS /path/to/chaste-install NO_DEFAULT_PATH)

# Here we add extra arguments to force chaste_codegen to use this extra argument
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
