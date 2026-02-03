#
# We need to treat separately the files in this folder to avoid duplication
# of symbols when building the tests
# If we do not do this and include these sources in the ${_this_exe} library,
# CUDA compilation of the tests will fail
#
set(MFIX_UDFS_STUBS
   ${CMAKE_CURRENT_LIST_DIR}/mfix_usr_drag.cpp
   ${CMAKE_CURRENT_LIST_DIR}/mfix_pc_usr.cpp
   CACHE INTERNAL "MFIX default UDFS")

target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_usr_reactions_rates_K.H
   )
