#
# Define a general build target to build all tests
#
set (BUILD_TEST build_tests)
add_custom_target(${BUILD_TEST})

#
# Define a macro to setup a test
#
macro (setup_test _test _suffix _env _input_file )

    string(REPLACE "/" "_" _test_name ${_test} )

    if (NOT TEST ${_test_name}_${_suffix})

        set(_test_dir  ${CMAKE_BINARY_DIR}/tests/${_test})
        set(_src_dir   ${PROJECT_SOURCE_DIR}/tests/${_test})

        # Copy file to binary dir
        file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${_test} DESTINATION ${CMAKE_BINARY_DIR}/tests/ )

        set(_test_exe mfix_${_test})

        # Create target if needed
        if (NOT TARGET ${_test_exe})
            add_executable(${_test_exe} EXCLUDE_FROM_ALL ${MAIN} )

            # Find UDFs if any is present and add them to the build
            file( GLOB _udfs  ${_src_dir}/*.cpp )
            add_udfs_to_target( ${_test_exe} ${_udfs} )

            target_include_directories(${_test_exe} PUBLIC ${_test_dir} )
            target_include_directories(${_test_exe} PRIVATE ${PROJECT_BINARY_DIR} )

            set_target_properties( ${_test_exe}
               PROPERTIES
               RUNTIME_OUTPUT_DIRECTORY ${_test_dir} )

            # Special treatment for CUDA
            if (MFIX_CUDA)
               setup_target_for_cuda_compilation(${_test_exe})
            endif ()

            target_link_libraries( ${_test_exe} mfixcore )

            # Accumulate the build for each test into target "build_tests"
            add_dependencies( ${BUILD_TEST} ${_test_exe} )

            # Create test dir in build path if needed
            # Update test files if needed
            add_custom_command(TARGET ${_test_exe} PRE_BUILD
                COMMAND mkdir -p ${_test_dir}
                COMMAND cp -R ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/* .
                WORKING_DIRECTORY ${_test_dir}
                COMMENT "Copying files to test run directory")
        endif ()

        # Setup for running the a single test
        set(_test_run run_${_test_name}_${_suffix})

        # FEXTRACT AND FJOIN_PAR are defined as internal cache variables
        add_custom_target(${_test_run} DEPENDS ${_test_exe}
          COMMAND ${CMAKE_COMMAND} -E env ${_env};FEXTRACT=${FEXTRACT};FJOIN_PAR=${FJOIN_PAR};
            ${_test_dir}/runtests.sh ${_test_dir}/${_test_exe} "\"\"" ${_input_file}
            WORKING_DIRECTORY  ${_test_dir}
            )

        # Make sure that fcompare and fextract are always built before a test is run
        add_dependencies(${_test_run} fjoin_par AMReX::fextract)

        add_test(NAME ${_test_name}_${_suffix}
            COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${_test_run}
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

    endif ()

endmacro()



############################################################
#                                                          #
#                      SGS TESTS                           #
#                                                          #
############################################################
set( SGS_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z
   FLD04-x FLD04-y FLD04-z
   FLD06-x FLD06-y FLD06-z
   FLD07-x FLD07-y FLD07-z
   FLD08-x FLD08-y FLD08-z
   FLD09-x FLD09-y FLD09-z
   DEM01-x DEM01-y DEM01-z
   DEM03-x DEM03-y DEM03-z
   DEM04-x DEM04-y DEM04-z
   DEM05-x DEM05-y DEM05-z
   DEM06-x DEM06-y DEM06-z
   DEM07
   monitors
   )


foreach( THIS_TEST ${SGS_TESTS})

    setup_test(
        ${THIS_TEST}
        "SGS"
        "MFIX_MPI=0;MFIX_OMP=0;OMP_NUM_THREADS=1"
        "inputs_single")

endforeach( THIS_TEST ${SGS_TESTS} )


############################################################
#                                                          #
#                      MGS & MGP TESTS                     #
#                                                          #
############################################################
set( MGS_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z
   FLD04-x FLD04-y FLD04-z
   FLD06-x FLD06-y FLD06-z
   FLD07-x FLD07-y FLD07-z
   FLD08-x FLD08-y FLD08-z
   FLD09-x FLD09-y FLD09-z
   DEM06-x DEM06-y DEM06-z
   DEM07
   )

foreach( THIS_TEST ${MGS_TESTS})

    setup_test(
        ${THIS_TEST}
        "MGS"
        "MFIX_MPI=0;MFIX_OMP=0;OMP_NUM_THREADS=1"
        "inputs_multiple")

    if (MFIX_MPI)
        setup_test(
            ${THIS_TEST}
            "MGP"
            "MFIX_MPI=1;MFIX_OMP=0;OMP_NUM_THREADS=1"
            "inputs_multiple")
    endif ()

endforeach( THIS_TEST ${MGS_TESTS} )


############################################################
#                                                          #
#                      TGS & TGSO TESTS                    #
#                                                          #
############################################################
set( TGS_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z
   FLD04-x FLD04-y FLD04-z
   FLD06-x FLD06-y FLD06-z
   FLD07-x FLD07-y FLD07-z
   FLD08-x FLD08-y FLD08-z
   FLD09-x FLD09-y FLD09-z
   DEM06-x DEM06-y DEM06-z
   DEM07
   )


foreach( THIS_TEST ${TGS_TESTS})

    setup_test(
        ${THIS_TEST}
        "TGS"
        "MFIX_MPI=0;MFIX_OMP=0;OMP_NUM_THREADS=1"
        "inputs_tiled")

    if (MFIX_OMP)
        setup_test(
            ${THIS_TEST}
            "TGSO"
            "MFIX_MPI=0;MFIX_OMP=1;OMP_NUM_THREADS=4"
            "inputs_tiled")
    endif ()

endforeach( THIS_TEST ${TGS_TESTS} )


############################################################
#                                                          #
#                      TGM TESTS                           #
#                                                          #
############################################################
if ( MFIX_OMP AND MFIX_MPI )

    set( TGM_TESTS
        DEM06-x DEM06-y DEM06-z
        )

    foreach( THIS_TEST ${TGM_TESTS})

        setup_test(
            ${THIS_TEST}
            "TGPO"
            "MFIX_MPI=1;MFIX_OMP=1;OMP_NUM_THREADS=2"
            "inputs_multiple_tiled")

    endforeach( )

endif()
