target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_set_mac_velocity_bcs.cpp
   mfix_compute_convective_term.cpp
   mfix_correct_small_cells.cpp
   mfix_compute_MAC_projected_velocities.cpp
   mfix_compute_forces.cpp
   )
