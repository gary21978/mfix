################################################################################
# CSG Parsing Tests
################################################################################

set(unit_tests_csg_SOURCES)
list(APPEND unit_tests_csg_SOURCES
  levelset/boolean.t.cpp
  levelset/empty.t.cpp
  levelset/extrude.t.cpp
  levelset/internal_flow.t.cpp
  levelset/primitives.t.cpp
  levelset/transform.t.cpp
  parser/boolean.t.cpp
  parser/empty.t.cpp
  parser/extrude.t.cpp
  parser/nest.cpp
  parser/other.t.cpp
  parser/primitives.t.cpp
  parser/transform.t.cpp
  )

if (CSG_CGAL_ENABLED)
   list(APPEND unit_tests_csg_SOURCES
     levelset/hull.t.cpp
     parser/hull.t.cpp
   )
endif()

add_executable(unit_tests_csg EXCLUDE_FROM_ALL
   ${unit_tests_csg_SOURCES}
  )

target_include_directories(unit_tests_csg
  PRIVATE ${CMAKE_SOURCE_DIR}/src
  )

find_package(Catch2 3 REQUIRED)
include(Catch)
target_link_libraries(unit_tests_csg
  Catch2::Catch2WithMain
  csg-eb
  )

if (CSG_CGAL_ENABLED)
   target_link_libraries(unit_tests_csg CGAL::CGAL)
endif()

catch_discover_tests(unit_tests_csg)
