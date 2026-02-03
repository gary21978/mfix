target_include_directories( ${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_stepdata.cpp
   mfix_leveldata.cpp
   mfix_fluid_update.cpp
   mfix_fluid_update_density.cpp
   mfix_fluid_update_species.cpp
   mfix_fluid_update_energy.cpp
   mfix_fluid_update_tracers.cpp
   mfix_fluid_update_velocity.cpp
   )
