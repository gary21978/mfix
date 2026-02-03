################################################################################
# CSG/EB Comparison Tests
################################################################################

add_executable(functional_tests
  benchmarks/05-cyl-fluidbed.cpp
  benchmarks/07-hopper.cpp
  benchmarks/08-cyclone.cpp
  benchmarks/02-settling.cpp
  benchmarks/tutorial-clr-prototype.cpp
  )

# Special treatment for CUDA in functional_tests
if(MFIX_CUDA)
   setup_target_for_cuda_compilation(functional_tests)
else()
   target_compile_features(functional_tests PRIVATE cxx_std_17)
endif()

get_filename_component(TUTORIAL_CLR_CSG "benchmarks/clr.csg" ABSOLUTE)

target_compile_definitions(functional_tests
                           PRIVATE TUTORIAL_CLR_CSG="${TUTORIAL_CLR_CSG}")

find_package(Git QUIET)
if(GIT_FOUND AND NOT EXISTS "${PROJECT_SOURCE_DIR}/subprojects/csg-eb/.git")
  message(FATAL_ERROR "CSG-EB not found.  Initialize the submodule by running:\n\n   git submodule update --init --recursive\n\n")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/subprojects/csg-eb)
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/subprojects/csg-eb)
find_package(Catch2 3 REQUIRED)

target_link_libraries(functional_tests
  AMReX::amrex
  Catch2::Catch2WithMain
  CsgEb::csg-eb
)
