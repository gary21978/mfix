target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_des_coll_K.H
   mfix_des_cond_K.H
   mfix_evolve_particles.cpp
   mfix_pc_init_collision_params.cpp
   mfix_pc_interactions_K.H
   )
