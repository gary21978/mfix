#include <mfix_eb.H>
#include <mfix_fix_inputs.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>

#include <csg.hpp>

using namespace amrex;

void MFIXEB::
make_eb_csg ( Vector<Geometry> const& a_geom )
{
  FixInputs fix("Nov. 2025");
  fix.swap<std::string>("mfix.geometry_filename", "csg.geometry_filename");

  bool is_internal_flow = true;
  Vector<Real> scaling_factor_vec(3, 1.0);
  Vector<Real> translation_vec(3, 0.0);

  ParmParse pp("csg");

  std::string csg_file;
  pp.query("geometry_filename", csg_file);

  if (csg_file.empty()) {

    reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Missing or invalid input: mfix.geometry_filename\n"
        << "Input required for mfix.geometry = csg";

  } else {
    Print() << "CSG geometry file: " << csg_file << '\n';
  }

  pp.query("internal_flow", is_internal_flow);

  if(pp.queryarr("scaling_factor", scaling_factor_vec, 0, 3)) {
    Print() << "WARNING: The implicit function magnitudes will not be scaled\n";
  }

  pp.queryarr("translation", translation_vec, 0, 3);

  Array<Real,3> scaling_factor = { scaling_factor_vec[0]
                                 , scaling_factor_vec[1]
                                 , scaling_factor_vec[2] };

  Array<Real,3> translation = { translation_vec[0]
                              , translation_vec[1]
                              , translation_vec[2] };

  if (is_internal_flow) { Print() << "\n Building geometry for internal flow\n"; }
  else { Print() << "\n Building geometry for external flow\n"; }

  auto csg_if = csg::get_csgif(csg_file, is_internal_flow);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(csg_if, "Unable to create CsgIF from geometry file");

  auto final_csg_if = EB2::translate(EB2::scale(*csg_if, scaling_factor), translation);

  auto gshop = EB2::makeShop(final_csg_if);
  build_levels(a_geom, gshop);
}
