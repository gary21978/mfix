target_include_directories(${_test_exe}
   PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_test_exe}
   PRIVATE
   chem_test.cpp
   chem_test_prob.cpp
   chem_test_check.cpp
   chem_test_check_eulerian01.cpp
   )

