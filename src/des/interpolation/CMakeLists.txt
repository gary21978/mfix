target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_interp_K.H
   mfix_eb_interp_K.H
   mfix_eb_interp_shepard_K.H
   mfix_eb_newton_solver_K.H
   mfix_eb_node_connectivity_K.H
   )
