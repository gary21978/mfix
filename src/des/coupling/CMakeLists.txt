target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_calc_transfer_coeffs.cpp
   mfix_coupling.cpp
   mfix_txfr_to_particle.cpp
   )
