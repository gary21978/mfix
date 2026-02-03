target_include_directories(${_this_exe}
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_eb.cpp
   mfix_eb_box.cpp
   mfix_eb_cylinder.cpp
   mfix_eb_hopper.cpp
   mfix_eb_regular.cpp
   mfix_eb_generic.cpp
   mfix_eb_stl.cpp
   mfix_eb_if.H
   mfix_write_eb_surface.cpp
   )

if ( MFIX_CSG )
   target_sources(${_this_exe} PRIVATE mfix_eb_csg.cpp)
   target_compile_definitions(${_this_exe} PUBLIC CSG_EB)
endif()
