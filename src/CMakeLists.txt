set(_this_exe "mfixcore")
add_library(mfixcore OBJECT)
#
# Add sources and header from current dir and all sub dirs
#

# Uncategorized
target_include_directories(mfixcore PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

# IMPORTANT: DO NOT ADD main.cpp HERE!
target_sources(mfixcore PRIVATE
   mfix.cpp
   #  main.cpp <---- DO NOT ADD THIS HERE
   mfix.H
   )

# I/0
add_subdirectory(io)
add_subdirectory(io/reports)

# Chemistry
add_subdirectory(chem)
add_subdirectory(chem/stiff_solvers)

# Convection
add_subdirectory(convection)

# DES
add_subdirectory(des)
add_subdirectory(des/dem)
add_subdirectory(des/pic)
add_subdirectory(des/deposition)
add_subdirectory(des/generator)
add_subdirectory(des/interpolation)
add_subdirectory(des/coupling)

# Diffusion
add_subdirectory(diffusion)

# Fluid
add_subdirectory(fluid)

# Regrid / FillPatch
add_subdirectory(regrid)

# Projection
add_subdirectory(projection)

# Setup
add_subdirectory(setup)

# Utility functions
add_subdirectory(util)

# Problem classes
add_subdirectory(prob)
add_subdirectory(prob/bc)
add_subdirectory(prob/distributions)

# User hooks (UDFs)
add_subdirectory(usr)

# Embedded boundary functions
add_subdirectory(eb)

# Time evolution functions
add_subdirectory(timestepping)

#
# Link to imported AMReX
#
target_link_libraries( mfixcore PUBLIC AMReX::amrex )

#
# Link to AMReX-Hydro
# Must include in sources too to enable propagation
# This is due to the fact that we are dealing with
# OBJECT libraries
#
target_sources(mfixcore PUBLIC $<TARGET_OBJECTS:amrex_hydro>) # Propagates objs
target_link_libraries( mfixcore PUBLIC amrex_hydro ) # Propagates includes

#
# Require C++14 or newer and no compiler extensions.
#
set_target_properties(mfixcore PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(mfixcore PUBLIC cxx_std_14)

#
# If MFIX_CUDA, set CUDA-specific properties and
# make C++ files be compiled as CUDA sources
#
if (MFIX_CUDA)
   setup_target_for_cuda_compilation(mfixcore)
endif ()

#
# Set the default flags
#
if ( NOT CMAKE_CXX_FLAGS )
   target_link_libraries(mfixcore PUBLIC AMReX::Flags_CXX)
endif ()

if (MFIX_FPE)
   target_link_libraries(mfixcore PUBLIC AMReX::Flags_FPE)
endif ()

#
# Add AMReX "build info" capabilities
#
include(AMReXBuildInfo)
generate_buildinfo(mfixcore ${CMAKE_SOURCE_DIR} REQUIRED)
