target_include_directories(${_this_exe} PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe} PRIVATE
   mfix_apply_nodal_projection.cpp
   mfix_proj_redistribution.cpp
   mfix_constraints.cpp
   mfix_mlmg_options.cpp
   )
