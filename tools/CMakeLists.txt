add_executable( fjoin_par EXCLUDE_FROM_ALL C++/fjoin_par.cpp )
target_link_libraries( fjoin_par AMReX::amrex )

if (MFIX_CUDA)
   setup_target_for_cuda_compilation(fjoin_par)
endif ()

# Define full path to fjoin_par
set(FJOIN_PAR ${CMAKE_CURRENT_BINARY_DIR}/fjoin_par CACHE INTERNAL "")



add_executable( pic2dem EXCLUDE_FROM_ALL PIC2DEM/main.cpp  )
target_link_libraries( pic2dem mfixcore )

add_subdirectory(PIC2DEM)

set_target_properties(pic2dem PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(pic2dem PUBLIC cxx_std_14)

if (MFIX_CUDA)
  setup_target_for_cuda_compilation(pic2dem)
endif ()

# Define full path to PIC2DEM
set(PIC2DEM ${CMAKE_CURRENT_BINARY_DIR}/PIC2DEM CACHE INTERNAL "")
