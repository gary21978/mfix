###############################################

# Here we define the default config options   #
# that can be overwritten by the user         #

###############################################
include_guard(GLOBAL)

include(CMakeDependentOption)

#
# Define a macro to print active options
#
macro (print_option _var)
   if (${_var})
      message( STATUS "   ${_var}")
   endif ()
endmacro ()


#
# Dimensionality of the build  ===============================================
#
set(AMReX_SPACEDIM_VALUES 1 2 3)
set(AMReX_SPACEDIM 3 CACHE STRING "Dimension of AMReX build: <1,2,3>")
list(REMOVE_DUPLICATES AMReX_SPACEDIM)
foreach(D IN LISTS AMReX_SPACEDIM)
    if(NOT D IN_LIST AMReX_SPACEDIM_VALUES)
       message(FATAL_ERROR "AMReX_SPACEDIM=${D} is not allowed."
          " Must be one of ${AMReX_SPACEDIM_VALUES}")
    endif()
endforeach()
message( STATUS "Building AMReX with AMReX_SPACEDIM = ${AMReX_SPACEDIM}")

#
# Populate the cache and check the value of the user-definable options
#
message(STATUS "Configuring AMReX with the following options enabled: ")

#
# This is the option to enable/disable xSDK mode
#
# To handle both XSDK options and corresponding plain AMReX options,
# we make use of policy CMP0077 introduced as a default in CMake 3.13
# Under policy CMP0077, normal variables prevent option()
# to set internal variables of the same name.
# Example: if XSDK mode is ON and XSDK_ENABLE_Fortran=ON, a normal
# variable AMReX_FORTRAN will be created and set to ON.
# This will stop the subsequent option( AMReX_FORTRAN "Enable Fortran language" ON )
# from being executed and no entry AMReX_FORTRAN will be created in the cache
#
option( USE_XSDK_DEFAULTS "Enable xSDK mode"  OFF )
print_option( USE_XSDK_DEFAULTS )


#
# Option to control the type of library: static vs shared  ===================
#
# Defaults:
#   CMake: static
#   xSDK: shared
# Precedence of user options:
#   AMReX_BUILD_SHARED_LIBS > BUILD_SHARED_LIBS > USE_XSDK_DEFAULTS
#   default:    unset              unset               OFF
#
get_property(SHARED_LIBS_SUPPORTED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
mark_as_advanced(BUILD_SHARED_LIBS)
if(DEFINED BUILD_SHARED_LIBS)
    set(AMReX_BUILD_SHARED_LIBS_DEFAULT ${BUILD_SHARED_LIBS})
elseif(USE_XSDK_DEFAULTS)
    set(AMReX_BUILD_SHARED_LIBS_DEFAULT ON)
else()
    set(AMReX_BUILD_SHARED_LIBS_DEFAULT OFF)
endif()
option(AMReX_BUILD_SHARED_LIBS "Build AMReX shared library"
    ${AMReX_BUILD_SHARED_LIBS_DEFAULT})
if(AMReX_BUILD_SHARED_LIBS AND NOT SHARED_LIBS_SUPPORTED)
    message(FATAL_ERROR "AMReX_BUILD_SHARED_LIBS requested but not supported by platform")
endif()
print_option( AMReX_BUILD_SHARED_LIBS )

#
# Option to control generation of install targets
#
option(AMReX_INSTALL "Generate Install Targets" ON)


#
# Option to control if Fortran must be enabled  ==============================
#
if ( USE_XSDK_DEFAULTS )
   option( XSDK_ENABLE_Fortran "Enable Fortran language" OFF )
   set( AMReX_FORTRAN ${XSDK_ENABLE_Fortran} )
   print_option(XSDK_ENABLE_Fortran)
endif ()

option( AMReX_FORTRAN "Enable Fortran language" OFF )
print_option( AMReX_FORTRAN )

#
# Option to control precision of the build  ==================================
#
if ( USE_XSDK_DEFAULTS )
   set( XSDK_PRECISION "DOUBLE" CACHE STRING "Precision:<SINGLE,DOUBLE>" )
   set( AMReX_PRECISION ${XSDK_PRECISION})
   set( AMReX_PARTICLES_PRECISION ${XSDK_PRECISION})
   message( STATUS "   XSDK_PRECISION = ${XSDK_PRECISION}")
endif ()


set(AMReX_PRECISION_VALUES SINGLE DOUBLE)
set(AMReX_PRECISION DOUBLE CACHE STRING "Precision of AMReX build: <SINGLE,DOUBLE>")
set_property(CACHE AMReX_PRECISION PROPERTY STRINGS ${AMReX_PRECISION_VALUES})
if(NOT AMReX_PRECISION IN_LIST AMReX_PRECISION_VALUES)
   message(FATAL_ERROR "AMReX_PRECISION=${AMReX_PRECISION} not supported."
      " Must be one of ${AMReX_PRECISION_VALUES}")
endif()
message( STATUS "   AMReX_PRECISION = ${AMReX_PRECISION}")


#
# GPU backends    =============================================================
#
set(AMReX_GPU_BACKEND_VALUES NONE CUDA)
set(AMReX_GPU_BACKEND NONE CACHE STRING "On-node, accelerated GPU backend: <NONE,CUDA>")
set_property(CACHE AMReX_GPU_BACKEND PROPERTY STRINGS ${AMReX_GPU_BACKEND_VALUES})
if (NOT AMReX_GPU_BACKEND IN_LIST AMReX_GPU_BACKEND_VALUES)
   message(FATAL_ERROR "AMReX_GPU_BACKEND=${AMReX_GPU_BACKEND} is not allowed."
      " Must be one of ${AMReX_GPU_BACKEND_VALUES}")
endif ()

if (NOT AMReX_GPU_BACKEND STREQUAL NONE)
   message( STATUS "   AMReX_GPU_BACKEND = ${AMReX_GPU_BACKEND}")

   # We might set different default for different GPUs in the future.
   set(AMReX_GPU_MAX_THREADS_DEFAULT "256")
   set(AMReX_GPU_MAX_THREADS ${AMReX_GPU_MAX_THREADS_DEFAULT} CACHE STRING
       "Maximum number of GPU threads per block" )
   message( STATUS "   AMReX_GPU_MAX_THREADS = ${AMReX_GPU_MAX_THREADS}")
endif ()

# Legacy variables for internal use only
if (AMReX_GPU_BACKEND STREQUAL CUDA)
   set(AMReX_CUDA  ON )
else ()
   set(AMReX_CUDA  OFF)
endif ()
set(AMReX_SYCL OFF)
set(AMReX_HIP  OFF)

#
# GPU RDC support
#
# https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_SEPARABLE_COMPILATION.html
set(_GPU_RDC_default ON)
if(AMReX_CUDA AND DEFINED CMAKE_CUDA_SEPARABLE_COMPILATION)
    set(_GPU_RDC_default "${CMAKE_CUDA_SEPARABLE_COMPILATION}")
endif()
cmake_dependent_option( AMReX_GPU_RDC "Enable Relocatable Device Code" ${_GPU_RDC_default}
   "AMReX_CUDA" OFF)
unset(_GPU_RDC_default)
print_option(AMReX_GPU_RDC)

#
# Fast Math    ================================================================
#
option(AMReX_FASTMATH  "Enable fast-math optimizations" OFF)
print_option(AMReX_FASTMATH)

#
# Parallel backends    ========================================================
#
option( AMReX_MPI  "Enable MPI"  ON )
print_option( AMReX_MPI )

cmake_dependent_option( AMReX_MPI_THREAD_MULTIPLE
   "whether to initialize MPI so that multiple threads can make MPI calls at the same time"  OFF
   "AMReX_MPI" OFF)
print_option( AMReX_MPI_THREAD_MULTIPLE )

option( AMReX_SIMD  "Enable SIMD Primitives" OFF)
print_option( AMReX_SIMD )

option( AMReX_OMP  "Enable OpenMP" OFF)
print_option( AMReX_OMP )


#
# AMReX components selection  ================================================
#
option( AMReX_AMRLEVEL  "Build AmrLevel class" ON )
print_option( AMReX_AMRLEVEL )

option( AMReX_EB "Build with Embedded Boundary support" OFF)
print_option(AMReX_EB)

cmake_dependent_option( AMReX_FORTRAN_INTERFACES "Build Fortran API" OFF
   "AMReX_FORTRAN" OFF )
print_option(AMReX_FORTRAN_INTERFACES)

option( AMReX_LINEAR_SOLVERS  "Build AMReX Linear solvers" ON )
print_option( AMReX_LINEAR_SOLVERS )

cmake_dependent_option( AMReX_LINEAR_SOLVERS_INCFLO
    "Build AMReX Linear solvers useful for incompressible flow codes" ON
    "AMReX_LINEAR_SOLVERS" OFF)
print_option( AMReX_LINEAR_SOLVERS_INCFLO )

cmake_dependent_option( AMReX_LINEAR_SOLVERS_EM
    "Build AMReX Linear solvers useful for electromagnetic codes" ON
    "AMReX_LINEAR_SOLVERS" OFF)
print_option( AMReX_LINEAR_SOLVERS_EM )

option( AMReX_AMRDATA "Build data services" OFF )
print_option( AMReX_AMRDATA )

option( AMReX_PARTICLES "Build particle classes" ON)
print_option( AMReX_PARTICLES )

if (AMReX_PARTICLES)
   set(AMReX_PARTICLES_PRECISION_VALUES SINGLE DOUBLE)
   set(AMReX_PARTICLES_PRECISION ${AMReX_PRECISION}
      CACHE STRING "Precision of reals in particle classes: <SINGLE,DOUBLE>")
   set_property(CACHE AMReX_PARTICLES_PRECISION PROPERTY STRINGS ${AMReX_PARTICLES_PRECISION_VALUES})
   if(NOT AMReX_PARTICLES_PRECISION IN_LIST AMReX_PARTICLES_PRECISION_VALUES)
      message(FATAL_ERROR "AMReX_PARTICLES_PRECISION=${AMReX_PRECISION} not supported."
         " Must be one of ${AMReX_PARTICLES_PRECISION_VALUES}")
   endif()
   message( STATUS "   AMReX_PARTICLES_PRECISION = ${AMReX_PARTICLES_PRECISION}")
endif ()


#
# External packages  =========================================================
#

# sensei
option( AMReX_SENSEI "Enable SENSEI in situ infrastructure" OFF )
print_option( AMReX_SENSEI )

cmake_dependent_option( AMReX_NO_SENSEI_AMR_INST
   "Disables the SENSEI instrumentation in amrex::Amr" FALSE
   "AMReX_SENSEI" FALSE )
print_option( AMReX_NO_SENSEI_AMR_INST )

# Conduit (requires CONDUIT_DIR)
option( AMReX_CONDUIT "Enable Conduit support" OFF )
print_option( AMReX_CONDUIT )

# Catalyst
cmake_dependent_option( AMReX_CATALYST "Enable Catalyst support" OFF 
   "AMReX_CONDUIT" OFF )
print_option( AMReX_CATALYST )

# Ascent
cmake_dependent_option( AMReX_ASCENT "Enable Ascent support" OFF
   "AMReX_CONDUIT" OFF )
print_option( AMReX_ASCENT )

# Hypre
cmake_dependent_option(AMReX_HYPRE "Enable Hypre interfaces" OFF
   "AMReX_LINEAR_SOLVERS" OFF)
print_option(AMReX_HYPRE)

# PETSc
cmake_dependent_option(AMReX_PETSC "Enable PETSc interfaces" OFF
   "AMReX_LINEAR_SOLVERS" OFF )
print_option(AMReX_PETSC)

# HDF5
option(AMReX_HDF5 "Enable HDF5-based I/O" OFF)
print_option(AMReX_HDF5)

cmake_dependent_option(AMReX_HDF5_ASYNC "Enable asynchronous writes in the HDF5-based IO" OFF
   "AMReX_HDF5" OFF )
print_option(AMReX_HDF5_ASYNC)

if (AMReX_HDF5_ASYNC)
   message(FATAL_ERROR "\nAMReX_HDF5_ASYNC not yet supported\n")
endif ()

cmake_dependent_option(AMReX_HDF5_ZFP "Enable ZFP compression in HDF5-based IO" OFF
   "AMReX_HDF5" OFF )
print_option(AMReX_HDF5_ZFP)

# SUNDIALS
option( AMReX_SUNDIALS "Enable SUNDIALS interfaces" OFF )
print_option( AMReX_SUNDIALS )


#
# Miscellanoues options  =====================================================
#
option( AMReX_PIC "Build position-independent code" OFF)
print_option( AMReX_PIC )

option( AMReX_IPO "Enable interprocedural optimization (IPO/LTO)" OFF)
print_option( AMReX_IPO )

option(AMReX_FPE "Enable Floating Point Exceptions checks" OFF)
print_option( AMReX_FPE )

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   option(AMReX_COMPILER_DEFAULT_INLINE "Use compiler default inline behavior" OFF)
   set(AMReX_INLINE_LIMIT 43210 CACHE STRING "Inline limit")
   if (NOT AMReX_COMPILER_DEFAULT_INLINE)
      if (AMReX_INLINE_LIMIT LESS 0)
         message(FATAL_ERROR "AMReX_INLINE_LIMIT, if set, must be non-negative")
      endif()
      message(STATUS "   AMReX_INLINE_LIMIT = ${AMReX_INLINE_LIMIT}")
   endif ()
else ()
   set(AMReX_COMPILER_DEFAULT_INLINE ON)
endif ()


if ( "${CMAKE_BUILD_TYPE}" MATCHES "Debug" )
   option( AMReX_ASSERTIONS "Enable assertions" ON)
else ()
   option( AMReX_ASSERTIONS "Enable assertions" OFF)
endif ()

print_option( AMReX_ASSERTIONS )

option( AMReX_FLATTEN_FOR "Enable flattening of ParallelFor and other similar functions" OFF)
print_option( AMReX_FLATTEN_FOR )

option(AMReX_BOUND_CHECK  "Enable bound checking in Array4 class" OFF)
print_option( AMReX_BOUND_CHECK )

if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    set(AMReX_EXPORT_DYNAMIC_DEFAULT ON)
else()
    set(AMReX_EXPORT_DYNAMIC_DEFAULT OFF)
endif()
option( AMReX_EXPORT_DYNAMIC "Enable Backtrace for macOS/Darwin" ${AMReX_EXPORT_DYNAMIC_DEFAULT})
print_option( AMReX_EXPORT_DYNAMIC )

#
# Profiling options  =========================================================
#
option( AMReX_BASE_PROFILE "Enable basic profiling" OFF )
print_option( AMReX_BASE_PROFILE )

cmake_dependent_option( AMReX_TINY_PROFILE "Enable tiny profiling" OFF
   "NOT AMReX_BASE_PROFILE" OFF)
print_option( AMReX_TINY_PROFILE )

cmake_dependent_option( AMReX_TRACE_PROFILE "Enable trace-profiling" OFF
   "AMReX_BASE_PROFILE" OFF)
print_option( AMReX_TRACE_PROFILE )

option( AMReX_MEM_PROFILE   "Enable memory profiling" OFF )
print_option( AMReX_MEM_PROFILE )

cmake_dependent_option( AMReX_COMM_PROFILE  "Enable communicator-profiling" OFF
   "AMReX_BASE_PROFILE" OFF)
print_option( AMReX_COMM_PROFILE )

cmake_dependent_option(AMReX_PROFPARSER "Enable profile parser" OFF
   "AMReX_BASE_PROFILE;AMReX_TRACE_PROFILE;AMReX_AMRDATA" OFF)
print_option( AMReX_PROFPARSER )


set(AMReX_TP_PROFILE_VALUES IGNORE CRAYPAT FORGE VTUNE)
set(AMReX_TP_PROFILE IGNORE CACHE STRING "Third-party profiling options: <CRAYPAT,FORGE,VTUNE>")
set_property(CACHE AMReX_TP_PROFILE PROPERTY STRINGS ${AMReX_TP_PROFILE_VALUES})
if (NOT AMReX_TP_PROFILE IN_LIST AMReX_TP_PROFILE_VALUES)
   message(FATAL_ERROR "AMReX_TP_PROFILE (${AMReX_TP_PROFILE}) must be one of ${AMReX_TP_PROFILE_VALUES}")
endif()
if (AMReX_TP_PROFILE)
   message( STATUS "   AMReX_TP_PROFILE = ${AMReX_TP_PROFILE}")
endif ()

# Check profile options
if ( NOT ( CMAKE_C_COMPILER_ID STREQUAL "Intel" ) AND
      ( AMReX_TP_PROFILE STREQUAL "VTUNE") )
   message( FATAL_ERROR "VTUNE cannot be used with ${CMAKE_C_COMPILER_ID} compiler" )
endif ()

if (  ( ( AMReX_TP_PROFILE STREQUAL "CRAYPAT" ) OR
        ( AMReX_TP_PROFILE STREQUAL "FORGE"   ) OR
        ( AMReX_TP_PROFILE STREQUAL "VTUNE"   )   ) AND
     (AMReX_BASE_PROFILE OR AMReX_TINY_PROFILE) )
   message(WARNING "This configuration should only be used to profile BL_PROFILE!")
endif()


#
# Extra options  =========================================================
#
option(AMReX_DIFFERENT_COMPILER
   "Allow an application to use a different compiler than the one used to build AMReX" OFF)
print_option(AMReX_DIFFERENT_COMPILER)

if ( NOT (CMAKE_SYSTEM_NAME STREQUAL "Linux") )
   option(AMReX_PROBINIT "Enable support for probin file" OFF)
else ()
   cmake_dependent_option(AMReX_PROBINIT "Enable support for probin file" ON
       "AMReX_AMRLEVEL" OFF)
endif ()
print_option(AMReX_PROBINIT)


#
# Static code analysis  ===============================================
#
option(AMReX_CLANG_TIDY "Enable clang-tidy analysis" OFF)
print_option(AMReX_CLANG_TIDY)
cmake_dependent_option(AMReX_CLANG_TIDY_WERROR "Treat clang-tidy warnings as errors" OFF
   "AMReX_CLANG_TIDY" OFF)
print_option(AMReX_CLANG_TIDY_WERROR)

#
# Tests  =============================================================
#
option(AMReX_ENABLE_TESTS "Enable CTest suite for AMReX" NO)
print_option(AMReX_ENABLE_TESTS)
set(AMReX_TEST_TYPE_VALUES "All;Small")
set(AMReX_TEST_TYPE All CACHE STRING "Type of AMReX Tests: <All,Small>")
set_property(CACHE AMReX_TEST_TYPE PROPERTY STRINGS ${AMReX_TEST_TYPE_VALUES})
if (NOT AMReX_TEST_TYPE IN_LIST AMReX_TEST_TYPE_VALUES)
   message(FATAL_ERROR "AMReX_TEST_TYPE=${AMReX_TEST_TYPE} is not allowed."
                       " Must be one of ${AMReX_TEST_TYPE_VALUES}.")
endif()
if (AMReX_ENABLE_TESTS)
   message(STATUS "   AMReX_TEST_TYPE = ${AMReX_TEST_TYPE}")
endif()

