target_include_directories( ${_this_exe} PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   integrator.H
   linear_algebra.H
   fe_integrator.H
   be_integrator.H
   vode_integrator.H
   )
