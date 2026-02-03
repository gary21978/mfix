target_include_directories(${_this_exe} PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe} PRIVATE
   mfix_init_fluid.cpp
   mfix_init.cpp
   mfix_init_redistribution.cpp
   mfix_set_gp0.cpp
   mfix_restart.cpp
   )
