target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_pc.cpp
   mfix_pc_remove_out_of_range.cpp
   mfix_pc_replicate.cpp
   mfix_pc_inflow.cpp
   mfix_pc_impose_mean.cpp
   )
