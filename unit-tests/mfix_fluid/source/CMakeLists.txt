target_include_directories(${_test_exe}
   PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_test_exe}
   PRIVATE
   test_fluid_incompressible.cpp
   test_fluid_molecular_viscosity.cpp
   test_fluid_suspension_viscosity.cpp
   test_fluid_eddy_viscosity.cpp
   test_fluid_porous_media.cpp
   )
