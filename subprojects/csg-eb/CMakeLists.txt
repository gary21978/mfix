cmake_minimum_required(VERSION 3.14)

project(CSG-EB
  DESCRIPTION  "Parser for CSG files to define Embedded Boundaries"
  HOMEPAGE_URL "https://mfix.netl.doe.gov/gitlab/exa/csg-eb"
  LANGUAGES    CXX
  )

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CSG_CGAL_ENABLED   "Enable CGAL"   ON)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set( CMAKE_CXX_COMPILER_LAUNCHER ccache )
endif()

find_package(pegtl REQUIRED)

if (CSG_CGAL_ENABLED)
   find_package(CGAL REQUIRED)
endif()

add_subdirectory(src/csg)

target_include_directories(csg-eb
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
  )

target_compile_features(csg-eb PRIVATE cxx_std_17)


file(GLOB headers "include/*.hpp")
install(FILES ${headers}
  DESTINATION include
  )

install(TARGETS csg-eb
    EXPORT CsgEbConfig
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
export(TARGETS csg-eb csg-eb-pi
    NAMESPACE CsgEb::
    FILE "${CMAKE_CURRENT_BINARY_DIR}/CsgEbConfig.cmake"
)
install(EXPORT CsgEbConfig
    DESTINATION "${CMAKE_INSTALL_PREFIX}"
    NAMESPACE CsgEb::
)


if (CSG_CGAL_ENABLED)
   target_compile_definitions(csg-eb PUBLIC USE_CGAL)
endif()

set(CSGEB_GIT_HASH "")
execute_process(
   COMMAND git describe --abbrev=12 --dirty --always --tags
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
   OUTPUT_VARIABLE CSGEB_GIT_HASH
   ERROR_VARIABLE  _tmperr
   OUTPUT_STRIP_TRAILING_WHITESPACE
   )
if(_tmperr)
   message(WARNING "Failing to retrieve Git commit from repo: ${PROJECT_SOURCE_DIR}")
endif()
set_target_properties(csg-eb PROPERTIES CSGEB_GIT_HASH "${CSGEB_GIT_HASH}"
   EXPORT_PROPERTIES CSGEB_GIT_HASH)

# to be used by gnumake
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CsgEbConfig.mk 
   "CSGEB_GIT_HASH := ${CSGEB_GIT_HASH}")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CsgEbConfig.mk
    DESTINATION "${CMAKE_INSTALL_PREFIX}"
)
