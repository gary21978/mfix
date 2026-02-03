target_include_directories( ${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_apply_predictor.cpp
   mfix_apply_corrector.cpp
   mfix_compute_dt.cpp
   mfix_deposit_to_grid.cpp
   mfix_evolve.cpp
   mfix_evolve_Godunov.cpp
   mfix_evolve_MOL.cpp
   mfix_evolve_fluid.cpp
   mfix_evolve_solids.cpp
   mfix_evolve_steady_state.cpp
   mfix_fluid_implicit_updates.cpp
   mfix_fluid_startup.cpp
   )
