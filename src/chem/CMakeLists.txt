target_include_directories(${_this_exe} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_actual_rhs.H
   mfix_calc_particles_within_cells.cpp
   mfix_chem_data_structs.H
   mfix_chem_memory_handler.H
   mfix_chem_memory_structs.H
   mfix_cmi.cpp
   mfix_cmi.H
   mfix_integrator.H
   mfix_material_interface.H
   )
