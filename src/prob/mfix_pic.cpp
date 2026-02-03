#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <mfix_reporter.H>
#include <mfix_pic.H>

using namespace amrex;


void
MFIXPIC::Initialize ()
{
  amrex::ParmParse pp("pic");

  // Names of the solids used to build input regions.
  amrex::Vector<std::string> names;
  pp.queryarr("solve", names);

  m_solve = names.size();

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_solve >= 0, "pic.m_solve must be >= 0");

  for(int lc=0; lc < names.size(); ++lc){
    if (amrex::toLower(names[0]).compare("none") == 0 ||
      (names[0]).compare("0") == 0) m_solve = 0;
  }

  // You can't name a solids "None" or "0" -- you just can't
  if (m_solve == 0 && names.size() > 1) {
    amrex::Abort("Invalid input: One or more PIC solids defined"
                  "but, the m_solver is disabled!");
  }

  if (m_solve) {

    // Store the total number of solids
    m_NPHASE = names.size();

    // verbosity level
    pp.query("verbose", m_verbose);

    // pic.pressure_coefficient ..................................................
    if(!pp.query("pressure_coefficient", m_Ps)) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input pic.pressure_coefficient not found!\n"
        << "Please correct the input deck.";

    } else if ( m_Ps <= 0. ) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.pressure_coefficient must be > 0\n"
        << "Please correct the input deck.";
    }

    // pic.beta ..................................................................
    if (!pp.query("beta", m_beta)) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input pic.beta not found!\n"
        << "Please correct the input deck.";

    } else if ( m_beta < 2. && m_beta > 5.) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.beta must be in range [2.0, 5.0]\n"
        << "Please correct the input deck.";
    }

    // pic.close_pack ............................................................
    if (!pp.query("close_pack", m_ep_cp)) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input pic.close_pack not found!\n"
        << "Please correct the input deck.";

    } else if ( m_ep_cp <= 0. || m_ep_cp >= 1.) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.close_pack must be in range (0.0, 1.0)\n"
        << "Please correct the input deck.";
    }

    // pic.parcels_per_cell_at_pack ..............................................
    if (!pp.query("parcels_per_cell_at_pack", m_parcels_per_cell_at_pack)) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input pic.parcels_per_cell_at_pack not found!\n"
        << "Please correct the input deck.";

    } else if ( m_parcels_per_cell_at_pack < 12 ) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.parcels_per_cell_at_pack must be > 12\n"
        << "Please correct the input deck.";
    }

    // pic.damping_factor ........................................................
    if (!pp.query("damping_factor", m_damping_factor)) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input pic.damping_factor not found!\n"
        << "Please correct the input deck.";

    } else if (m_damping_factor < 0. && m_damping_factor > 1.) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.damping_factor must be in range [0.0, 1.0]\n"
        << "Please correct the input deck.";
    }

    // pic.damping_factor_wall_normal ............................................
    if (!pp.query("damping_factor_wall_normal", m_damping_factor_wall_normal)) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input pic.damping_factor_wall_normal not found!\n"
        << "Please correct the input deck.";

    } else if (m_damping_factor_wall_normal < 0. || m_damping_factor_wall_normal > 1.) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.damping_factor_wall_normal must be in range [0.0, 1.0]\n"
        << "Please correct the input deck.";
    }

    // pic.damping_factor_wall_tangent ...........................................
    if (!pp.query("damping_factor_wall_tangent", m_damping_factor_wall_tangent)) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input pic.damping_factor_wall_tangent not found!\n"
        << "Please correct the input deck.";

    } else if (m_damping_factor_wall_tangent < 0. || m_damping_factor_wall_tangent > 1.) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.damping_factor_wall_tangent must be in range [0.0, 1.0]\n"
        << "Please correct the input deck.";
    }

    // pic.damping_factor_wall_tangent ...........................................
    pp.query("small_number", m_small_number);
    if (m_small_number <= 0.) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.small_number must be > 0\n"
        << "Please correct the input deck.";
    }

    // pic.velocity_reference_frame ..............................................
    pp.query("velocity_reference_frame", m_vel_ref_frame);
    if (m_vel_ref_frame < 0. || m_vel_ref_frame > 1.) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.velocity_reference_frame must be in range [0.0, 1.0]\n"
        << "Please correct the input deck.";
    }

    // pic.advance_vel_p .........................................................
    pp.query("advance_vel_p", m_advance_vel_p);
    if (m_advance_vel_p < 0.0 || m_advance_vel_p > 1.0) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.advance_vel_p must be in range [0.0, 1.0]\n"
        << "Please correct the input deck.";
    }

    // pic.max_iter ..............................................................
    pp.query("max_iter", m_max_iter);
    if (m_max_iter < 1) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.max_iter must be >= 1\n"
        << "Please correct the input deck.";
    }

    // pic.initial_step_type .....................................................
    std::string l_initial_step = "nth_eps";
    pp.queryAdd("initial_step_type", l_initial_step);
    l_initial_step = toLower(l_initial_step); // case insensitive

    if (l_initial_step == "nth_eps") {
      m_initial_step = InitialStepType::nth_eps;
    } else if (l_initial_step == "zero_eps") {
      m_initial_step = InitialStepType::zero_eps;
    } else if (l_initial_step == "taylor_approx") {
      m_initial_step = InitialStepType::taylor_approx;
    } else {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid value: pic.initial_step_type must be nth_eps, zero_eps or taylor_approx\n"
        << "Please correct the input deck.";
    }

  }
}
