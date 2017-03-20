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

# Here we just depend on core components (nothing application-specific).
find_package(Chaste COMPONENTS heart)

# Alternatively, to specify a Chaste installation directory use a line like that below.
# This is needed if your project is not contained in the projects folder within a Chaste source tree.
#find_package(Chaste COMPONENTS heart crypt PATHS /path/to/chaste-install NO_DEFAULT_PATH)

# Here we add extra arguments to force PyCML to use this extra argument (make Get and Set methods for
# all metadata annotated variables). 
set(Chaste_PYCML_EXTRA_ARGS "--expose-annotated-variables")

# ApPredict uses pthreads, so we specify that threads are required and add the necessary compile flag
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
if(CMAKE_USE_PTHREADS_INIT)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -pthread")
endif()

# Change the project name in the line below to match the folder this file is in,
# i.e. the name of your project.
chaste_do_project(ApPredict)
