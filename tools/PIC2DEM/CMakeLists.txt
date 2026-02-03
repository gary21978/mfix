target_include_directories( pic2dem
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( pic2dem
   PRIVATE
   restarter.H
   restarter.cpp
   mfix_depop_restarter.cpp
   )

 file( GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.cpp )
 add_udfs_to_target(pic2dem ${USR_OVERRIDES})
