target_include_directories(${_this_exe} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( ${_this_exe}
   PRIVATE
   mfix_avg_pc_parms.cpp
   mfix_bc.cpp
   mfix_dem.cpp
   mfix_fluid.cpp
   mfix_ic.cpp
   mfix_load_balance.cpp
   mfix_pic.cpp
   mfix_reactions.cpp
   mfix_regions.cpp
   mfix_solids.cpp
   mfix_species.cpp
   mfix_tag.cpp
   mfix_timer.cpp
   mfix_porous_media.cpp
   mfix_properties.cpp
   mfix_prop_molecular_viscosity.cpp
   mfix_prop_molecular_weight.cpp
   mfix_prop_specific_heat.cpp
   )
