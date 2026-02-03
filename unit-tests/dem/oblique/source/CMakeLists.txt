target_include_directories(${_test_exe}
   PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_test_exe}
   PRIVATE
   dem_test.cpp
   dem_test_eb.cpp
   dem_test_prob.cpp
   dem_test_check.cpp
   mfix_pc_usr.cpp
   )
