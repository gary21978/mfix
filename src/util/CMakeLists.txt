target_include_directories(${_this_exe} PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe} PRIVATE
   mfix_filcc.cpp
   mfix_mf_helpers.cpp
   mfix_calc_cell.cpp
   )
