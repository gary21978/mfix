cmake_minimum_required(VERSION 3.20)

project(MFIX-Exa
   # VERSION <version>
   DESCRIPTION  "A multiphase modeling tool for the exascale"
   HOMEPAGE_URL "https://mfix.netl.doe.gov/doc/mfix-exa/guide/latest/"
   LANGUAGES    C CXX
)

# Set the search path for cmake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/tools/CMake)

# Always include debug info so we get useful stack traces
# Don't use -O3
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g -DNDEBUG")

# Load utilities and MFIX options module
include( MFIX_Utils )

# Creates `compile_commands.json` in the build build directory
#  `compile_commands.json` contains compiler flags used by plugins like YouCompleteMe
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# Set default build type to Release
if ( NOT CMAKE_BUILD_TYPE )
   message(STATUS "Setting build type to Release as none was specified.")
   set( CMAKE_BUILD_TYPE Release )
else ()
   message(STATUS "Build type set by user to '${CMAKE_BUILD_TYPE}'.")
endif()

find_program(LLD_FOUND lld)
if(LLD_FOUND)
  add_link_options("-fuse-ld=lld")
endif()

#find_program(CCACHE_FOUND ccache)
#find_program(SCCACHE_FOUND sccache)
#if(SCCACHE_FOUND AND NOT SKIP_LAUNCHER)
#  set( CMAKE_CXX_COMPILER_LAUNCHER sccache )
#  set( CMAKE_CUDA_COMPILER_LAUNCHER sccache )
#elseif(CCACHE_FOUND AND NOT SKIP_LAUNCHER)
#  set( CMAKE_CXX_COMPILER_LAUNCHER ccache )
#  set( CMAKE_CUDA_COMPILER_LAUNCHER ccache )
#endif()

# We want to offer the user some core options, i.e. MPI, OpenMP, CUDA, HYPRE
# and floating-point exceptions, which require related AMReX components to
# be enabled in the installation.
# In STANDALONE mode, find_package will make sure these options are available
# IN SUPERBUILD mode, These options will supersede AMReX options with the same name
option( MFIX_OMP    "Enable OpenMP" NO )
option( MFIX_MPI    "Enable MPI"   YES )

option( MFIX_POLYDISPERSE "Enable polydisperse support" NO )
option( MFIX_MPMD "Enable multiple program multiple data. Requires MFIX_MPI enabled" NO )

include(CMakeDependentOption)
cmake_dependent_option( MFIX_MPI_THREAD_MULTIPLE
   "whether to initialize MPI so that multiple threads can make MPI calls at the same time"  OFF
   "MFIX_MPI" OFF)

option( MFIX_HYPRE  "Enable HYPRE"  NO )
option( MFIX_FPE    "Enable Floating Point Exceptions checks" NO )

#
# GPU backends    =============================================================
#
set(MFIX_GPU_BACKEND_VALUES NONE SYCL CUDA HIP)
set(MFIX_GPU_BACKEND NONE CACHE STRING "On-node, accelerated GPU backend: <NONE,SYCL,CUDA,HIP>")
set_property(CACHE MFIX_GPU_BACKEND PROPERTY STRINGS ${MFIX_GPU_BACKEND_VALUES})
if (NOT MFIX_GPU_BACKEND IN_LIST MFIX_GPU_BACKEND_VALUES)
   message(FATAL_ERROR "MFIX_GPU_BACKEND=${MFIX_GPU_BACKEND} is not allowed."
      " Must be one of ${MFIX_GPU_BACKEND_VALUES}")
endif ()

if (NOT MFIX_GPU_BACKEND STREQUAL NONE)
   message( STATUS "   MFIX_GPU_BACKEND = ${MFIX_GPU_BACKEND}")
endif ()

if (MFIX_GPU_BACKEND STREQUAL SYCL)
   set(MFIX_DPCPP ON  )
   set(MFIX_CUDA  OFF )
   set(MFIX_HIP   OFF )
elseif (MFIX_GPU_BACKEND STREQUAL CUDA)
   set(MFIX_DPCPP OFF )
   set(MFIX_CUDA  ON  )
   set(MFIX_HIP   OFF )
elseif (MFIX_GPU_BACKEND STREQUAL HIP)
   set(MFIX_DPCPP OFF )
   set(MFIX_CUDA  OFF )
   set(MFIX_HIP   ON  )
else ()
   set(MFIX_DPCPP OFF )
   set(MFIX_CUDA  OFF )
   set(MFIX_HIP   OFF )
endif ()

# CSG support
option( MFIX_CSG "Build with CSG support" OFF)

# Global switch to turn on/off building of tests
option( MFIX_TESTS "Build MFIX tests" ON)

# Switch to turn on/off building of CSG tests
option( MFIX_CSG_TESTS "Build MFIX CSG tests" OFF)

# VODE support
option( VODE_JACOBIAN_CACHING "Build VODE integrator with Jacobian caching support" OFF)

if (VODE_JACOBIAN_CACHING)
  add_compile_definitions(VODE_JACOBIAN_CACHING)
endif()

if (MFIX_CUDA)
   enable_language(CUDA)
endif ()

# Setup AMReX
include(MFIXSetupAMReX)

# Setup AMReX-Hydro
include(MFIXSetupAMReXHydro)

# Setup CSG-EB
include(MFIXSetupCsgeb)

get_git_info( )

# Set some cache variable for the whole project
set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp CACHE INTERNAL "Main source file")

#
# Core MFIX library
#
add_subdirectory(src)

# Print out config summary
print_mfix_configuration_summary(mfixcore)

add_executable(mfix ${MAIN})
target_link_libraries( mfix mfixcore )
target_include_directories( mfix PRIVATE ${PROJECT_BINARY_DIR})
target_include_directories( mfixcore PRIVATE ${PROJECT_BINARY_DIR})

# We should check what this line does
file( GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.cpp)
add_udfs_to_target(mfix ${USR_OVERRIDES})

if (MFIX_CUDA)
   setup_target_for_cuda_compilation(mfix)
endif ()

# Tools
add_subdirectory(tools)

# Tests
if (MFIX_TESTS)
   enable_testing()
   add_subdirectory( tests )
   add_subdirectory( unit-tests )
   if (MFIX_CSG AND MFIX_CSG_TESTS )
      add_subdirectory(src/eb/tests)
   endif()
endif ()

# Tags
include(MFIX_Tags)
add_tags_targets()

# Copy compile_commands.json from the build directory to the project root directory
if ( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
    execute_process( COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
endif ()

if(MFIX_CSG)
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_BINARY_DIR}/subprojects/csg-eb)
  list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR}/subprojects/csg-eb)
  target_link_libraries(mfixcore PUBLIC CsgEb::csg-eb)
endif()

if(MFIX_POLYDISPERSE)
   target_compile_definitions(mfixcore PUBLIC MFIX_POLYDISPERSE)
endif()

if(MFIX_MPMD)
   if(NOT MFIX_MPI)
      message(FATAL_ERROR "MFIX_MPMD=${MFIX_MPMD} is not allowed when MFIX_MPI=${MFIX_MPI}.")
   endif()
   target_compile_definitions(mfixcore PUBLIC MFIX_MPMD)
endif()
