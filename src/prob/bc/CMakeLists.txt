target_include_directories(${_this_exe} PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe} PRIVATE
   mfix_calc_bc_areas.cpp
   mfix_set_bcs.cpp
   mfix_set_bc_list.cpp
   mfix_set_bc_type.cpp
   mfix_set_bc_values.cpp
   mfix_set_eb_scalar_bcs.cpp
   mfix_set_eb_temperature_bcs.cpp
   mfix_set_eb_velocity_bcs.cpp
   mfix_set_gradp_bcs.cpp
   mfix_bc_fillpatch.cpp
   mfix_bc_fillcoarsepatch.cpp
   )
