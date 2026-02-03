target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_diffusion_op.cpp
   mfix_diffusion_op.H
   mfix_diffop_energy.cpp
   mfix_diffop_scalar.cpp
   mfix_diffop_species.cpp
   mfix_diffop_tracer.cpp
   mfix_diffop_velocity_K.H
   mfix_diffop_velocity.cpp
   )
