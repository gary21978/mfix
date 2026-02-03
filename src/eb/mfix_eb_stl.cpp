#include <mfix_eb.H>
#include <mfix_fix_inputs.H>

using namespace amrex;

void MFIXEB::
make_eb_stl ( Vector<Geometry> const& a_geom )
{
  FixInputs fix("Nov. 2025");
  fix.swap<std::string>("mfix.geometry_filename", "stl.geometry_filename");

  bool is_internal_flow = true;
  Real scaling_factor = 1.0;
  Vector<Real> translation_vec(3, 0.0);
  bool use_bvh = true;

  ParmParse pp("stl");

  std::string stl_file;
  pp.query("geometry_filename", stl_file);

  if (stl_file.empty()) {

    reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Missing or invalid input: stl.geometry_filename\n"
        << "Input required for mfix.geometry = stl";

  } else {

    Print() << "STL geometry file: " << stl_file << '\n';

  }

  pp.query("internal_flow", is_internal_flow);
  pp.query("scaling_factor", scaling_factor);
  pp.queryarr("translation", translation_vec, 0, 3);
  pp.query("use_bvh", use_bvh);

  if (is_internal_flow) { Print() << "\n Building geometry for internal flow\n"; }
  else { Print() << "\n Building geometry for external flow\n"; }

  build_levels(a_geom, stl_file, scaling_factor,
      {translation_vec[0], translation_vec[1], translation_vec[2]},
      is_internal_flow, use_bvh);
}
