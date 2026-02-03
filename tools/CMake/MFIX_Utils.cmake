#
# Find git info
#
macro( get_git_info ) # EXTRA ARGS: branch commit

  # Find branch
  execute_process(
    COMMAND git branch
    COMMAND grep \*
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE out
    ERROR_VARIABLE  err
    )

  if(err)
    message(WARNING "Failing to retrieve MFIX Git branch")
  else()
    string( REPLACE "*" "" out ${out} )
    string( STRIP ${out} out )
    message( STATUS "MFIX branch: ${out}" )
    if( ${ARGC} GREATER 0 ) # branch
      set( ${ARGV0} ${out} )
    endif()
  endif()

  get_hash(${PROJECT_SOURCE_DIR}/subprojects/AMReX-Hydro)
  set(HYDRO_GIT_HASH ${GIT_HASH})

  configure_file(
	  tools/CMake/build_info.H.in
	  ${PROJECT_BINARY_DIR}/build_info.H
	  @ONLY)

  unset(out)

endmacro()

function( get_hash repo_path )
  execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${repo_path}
    OUTPUT_VARIABLE out
    ERROR_VARIABLE  err
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(err)
    message(WARNING "Failing to retrieve Git commit from repo: ${repo_path}")
  endif()
  set(GIT_HASH "${out}" PARENT_SCOPE)
endfunction()

#
# Print Configuration Summary
#
function( print_mfix_configuration_summary mfix_libname )

  if(NOT TARGET ${mfix_libname})
    message(AUTHOR_WARNING "Target ${mfix_libname} is not defined.")
    return()
  endif()

  if(NOT TARGET AMReX::amrex)
    message(AUTHOR_WARNING "Target ${mfix_libname} is not defined.")
    return()
  endif()

  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPERCASE)


  # Get AMReX cmake functions
  include(AMReXGenexHelpers)
  include(AMReXTargetHelpers)

  get_target_properties_flattened(${mfix_libname} _includes _defines _flags _link_line)

  set(_lang CXX)
  set(_prop _includes _defines _flags _link_line )


  # Loop over all combinations of language and property and extract
  # what you need
  foreach( _p IN LISTS _prop )
     foreach( _l IN LISTS _lang )

        string(TOLOWER ${_l} _ll) # Lower case language name

        # _${_ll}${_p} is a variable named as _lang_property,
        # both lower case.
        set(_${_ll}${_p} "${${_p}}")
        eval_genex( _${_ll}${_p} ${_l} ${CMAKE_${_l}_COMPILER_ID}
           COMP_VERSION ${CMAKE_${_l}_COMPILER_VERSION}
           CONFIG       ${CMAKE_BUILD_TYPE}
           INTERFACE    BUILD)

        if (_${_ll}${_p})

           list(REMOVE_DUPLICATES _${_ll}${_p})

           if ("${_p}" STREQUAL "_defines")
              string(REPLACE ";" " -D" _${_ll}${_p} "-D${_${_ll}${_p}}")
           elseif ("${_p}" STREQUAL "_includes")
              string(REPLACE ";" " -I" _${_ll}${_p} "-I${_${_ll}${_p}}")
           else()
              string(REPLACE ";" " " _${_ll}${_p} "${_${_ll}${_p}}")
           endif ()

        endif ()

     endforeach()
  endforeach ()


  #
  # Config summary
  #
  message( STATUS "MFIX configuration summary: " )
  message( STATUS "   C++ defines           = ${_cxx_defines}" )
  message( STATUS "   C++ compiler          = ${CMAKE_CXX_COMPILER}" )
  message( STATUS "   C++ flags             = ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}} ${CMAKE_CXX_FLAGS} ${_cxx_flags}" )
  if (MFIX_CUDA)
     message( STATUS "   CUDA flags            = ${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}} ${CMAKE_CUDA_FLAGS}" )
  endif ()
  message( STATUS "   MFIX includes         = ${_cxx_includes}" )
  message( STATUS "   MFIX link line        = ${_cxx_link_line}" )

endfunction()


#
# Adds a list of UDFs to _target source list
# If none is provided, it adds MFIX default UDFS stubs
#
function (add_udfs_to_target _target)

   if (NOT TARGET ${_target})
      message(FATAL_ERROR "Target ${_target} does not exist")
   endif ()

   set(_udfs ${ARGN})

   if (_udfs)
      set(_udfs_stubs ${MFIX_UDFS_STUBS})
      foreach (_source IN LISTS _udfs)
         get_filename_component(_source_name ${_source} NAME)
         list(FILTER _udfs_stubs EXCLUDE REGEX ${_source_name})
      endforeach ()
      list(APPEND _udfs ${_udfs_stubs})
   else ()
      list(APPEND _udfs ${MFIX_UDFS_STUBS})
   endif ()

   target_sources(${_target} PUBLIC ${_udfs})

endfunction ()
