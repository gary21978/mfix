target_include_directories( ${_this_exe} PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe} PRIVATE
   mfix_io.cpp
   mfix_build_info.cpp
   mfix_catalyst.cpp
   mfix_chk.cpp
   mfix_plt.cpp
   mfix_ascent.cpp
   mfix_rw.cpp
   mfix_monitors.cpp
   mfix_compute_vorticity.cpp
   )

if (MFIX_MPMD)
   target_sources( ${_this_exe} PRIVATE
      mfix_mpmd.cpp)
endif ()
