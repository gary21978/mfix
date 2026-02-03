set(_Bdirs des
           des/generator
           des/dem
           prob
           prob/distributions
           util
   )

set(_this_exe ${_test_exe})

foreach(_dir ${_Bdirs})
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/${_dir} ${_test_dir}/${_dir} EXCLUDE_FROM_ALL)
endforeach( )

add_subdirectory(source)

