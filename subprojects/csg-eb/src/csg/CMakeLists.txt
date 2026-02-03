################################################################################
# CSG Parser
################################################################################

set(csg-eb_SOURCES)
list(APPEND csg-eb_SOURCES
  csg.cpp
  levelset_2d.cpp
  levelset_3d.cpp
  matrix_functions.cpp
  parser.cpp
  exception.cpp
  )

if (CSG_CGAL_ENABLED)
   list(APPEND csg-eb_SOURCES
      cgal_helper_polygon.cpp
      cgal_helper_polyhedron.cpp
   )
endif()

## Create a private interface library
add_library(csg-eb-pi INTERFACE)

target_link_libraries(csg-eb-pi INTERFACE
                      stdc++fs
                      taocpp::pegtl
                      )
if (CSG_CGAL_ENABLED)
   target_link_libraries(csg-eb-pi INTERFACE CGAL::CGAL)
endif()

# Create a static library
add_library(csg-eb ${csg-eb_SOURCES})
add_library(CsgEb::csg-eb ALIAS csg-eb)
target_link_libraries(csg-eb PRIVATE "$<BUILD_INTERFACE:csg-eb-pi>")

if (BUILD_TESTING)
   add_subdirectory(tests)
endif()
