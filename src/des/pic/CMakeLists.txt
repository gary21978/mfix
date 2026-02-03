target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_evolve_parcels.cpp
   mfix_calc_avg_solids_prop.cpp
   mfix_calc_solids_stress.cpp
   mfix_pic_iteration.cpp
   mfix_pc_advance_pic.cpp
   mfix_pc_deposit_pic.cpp
   mfix_pic_K.H
   )
