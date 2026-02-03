target_include_directories(${_test_exe}
   PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_test_exe}
   PRIVATE
   test_constant_dist.cpp
   test_custom_dist.cpp
   test_dist.cpp
   test_log_normal_dist.cpp
   test_normal_dist.cpp
   test_uniform_dist.cpp
   )
