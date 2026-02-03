target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(${_this_exe}
   PRIVATE
   mfix_deposition_op.cpp
   mfix_depop_bcs.cpp
   mfix_depop_txfr.cpp
   mfix_depop_volavg.cpp
   mfix_depop_redist.cpp
   )
