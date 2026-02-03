#
# Define a general build target to build all tests
#
add_custom_target(unit_tests)

macro (setup_unit_test _test_name _test_base _inputs_file)

    if (NOT TEST ${_test_name})

        set(_test_dir ${CMAKE_BINARY_DIR}/unit-tests/${_test_base})
        set(_src_dir  ${PROJECT_SOURCE_DIR}/unit-tests/${_test_base})

        set(_test_exe ${_test_name})

        # Create target if needed
        if (NOT TARGET ${_test_exe})

            add_executable(${_test_exe} EXCLUDE_FROM_ALL ${_src_dir}/source/main.cpp )

            set_target_properties( ${_test_exe}
               PROPERTIES
               RUNTIME_OUTPUT_DIRECTORY ${_test_dir} )

            target_include_directories(${_test_exe} PUBLIC ${_test_dir} )
            target_include_directories(${_test_exe} PRIVATE ${PROJECT_BINARY_DIR} )

            add_subdirectory(${_src_dir})

            # If MFIX_CUDA, set CUDA-specific properties and compile C++ files as CUDA sources
            if (MFIX_CUDA)
               setup_target_for_cuda_compilation(${_test_exe})

               find_library(CURAND curand PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
               find_library(CUSPARSE cusparse PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
               target_link_libraries(${_test_exe} PUBLIC AMReX::amrex ${CURAND})

            else ()

                target_link_libraries( ${_test_exe} AMReX::amrex )

            endif ()

            # Accumulate the build for each test into target "unit_tests"
            add_dependencies( unit_tests ${_test_exe} )

            # Set test specific compile definitions
            foreach(_arg ${ARGN})
                target_compile_definitions(${_test_exe} PRIVATE ${_arg})
            endforeach()

            if (NOT ${_inputs_file} MATCHES "NONE")
                file( COPY ${PROJECT_SOURCE_DIR}/unit-tests/${_test_base}/${_inputs_file} DESTINATION ${_test_dir} )
                add_test(NAME unit_test-${_test_name}
                         COMMAND $<TARGET_FILE:${_test_exe}> "${_test_dir}/${_inputs_file}")
            else()
                add_test(NAME unit_test-${_test_name}
                         COMMAND $<TARGET_FILE:${_test_exe}>)
            endif()

        endif ()

    endif ()

endmacro()

############################################################################################
#                                                                                          #
#                                  Define unit tests                                       #
#                                                                                          #
############################################################################################
# Dummy arguments:
# 1 (required) : a unique name for the test
# 2 (required) : directory
# 3 (required) : name of inputs file :: set to NONE if an inputs is not needed
# 4 (optional) : compiler definitions
setup_unit_test(distributions distributions NONE)
setup_unit_test(mfix_fluid_var mfix_fluid_var NONE)
setup_unit_test(mfix_fluid mfix_fluid NONE MFIX_NO_ABORT MFIX_SILENT_REPORTER)
setup_unit_test(chemistry.eulerian.01 chemistry inputs.eulerian.01 CHEM_TEST_EULERIAN01)
setup_unit_test(dem.oblique dem/oblique inputs PARTICLE_INPUT_FILE=${PROJECT_SOURCE_DIR}/unit-tests/dem/oblique/particle_input.dat)
