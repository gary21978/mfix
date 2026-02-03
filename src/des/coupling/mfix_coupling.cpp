#include <AMReX.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <mfix_reporter.H>
#include <mfix_coupling.H>
#include <mfix_fix_inputs.H>

using namespace amrex;

int
CouplingOp::
Initialize ( int const a_couple_energy,
             int const a_nspecies,
             int const a_is_mixture,
             int const a_solve_fluid_and_solids )
{
  m_interp_idxs = new DualGridAuxIndexes(a_couple_energy, a_nspecies, a_is_mixture);

  // If we're not solving for both fluid and solids phase, we don't do anything
  // in here and we return immediately
  if (!a_solve_fluid_and_solids)
    return 0;

  // Fix old way for specifying drag model
  FixInputs fix("Oct. 2024");
  fix.swap<std::string>("mfix.drag_type", "mfix.drag.model");
  fix.swap<std::string>("mfix.SyamOBrien.c1", "mfix.drag.model.SyamOBrien.c1");
  fix.swap<std::string>("mfix.SyamOBrien.d1", "mfix.drag.model.SyamOBrien.d1");
  fix.swap<std::string>("mfix.convection_type", "mfix.convection.model");

  ParmParse pp_drag("mfix.drag");

  { std::string drag_model_str = "none";
    pp_drag.query("model", drag_model_str);

    std::string mlower(toLower(drag_model_str));
    if (mlower == "wenyu") {

      m_drag_model.setWenYu();

    } else if (mlower == "gidaspow") {

      m_drag_model.setGidaspow();

    } else if (mlower == "bvk2") {

      m_drag_model.setBVK2();

    } else if (mlower == "syamobrien") {

      m_drag_model.setSyamOBrien();

      ParmParse ppSyamOBrien("mfix.drag.model.SyamOBrien");
      Real c1, d1;
      if (ppSyamOBrien.query("c1", c1) && ppSyamOBrien.query("d1", d1) ) {

        m_drag_model.setSyamOBrien_coeffs(c1, d1);

      } else {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "SyamOBrien drag model requires coefficients c1 and d1.\n"
          << "Specify the following entries in the inputs file:\n"
          << "   mfix.drag.model.SyamOBrien.c1 = amrex::Real\n"
          << "   mfix.drag.model.SyamOBrien.d1 = amrex::Real";
        return 1;
      }

    } else if (mlower == "userdrag") {

      m_drag_model.setUserDrag();

    } else {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Unknown drag model: " << drag_model_str << "\n"
        << "Please correct the input deck.";
      return 1;
    }
    reporter::Log(reporter::Status) << "Drag model: " << drag_model_str;
  }

  // Include viscous stress in drag force on particles. The 'equal and opposite'
  // force is applied to the fluid by subtracting eps*div(tau) from the fluid.
  // The result :: div(tau) - eps*div(tau) = epg*div(tau).
  { int drag_include_divtau(0);
    if (pp_drag.query("include_divtau", drag_include_divtau)) {
      set_include_divtau(drag_include_divtau);
    }
    if (include_divtau()) {
      reporter::Log(reporter::Status) << "Drag forces includes div(tau)";
    }
  }

  // Include virtual mass force in drag force on particles.
  { std::string virtual_mass_model_str = "none";
    if (pp_drag.query("virtual_mass", virtual_mass_model_str)) {

      std::string mlower(toLower(virtual_mass_model_str));
      if (mlower == "null") {

            m_virtual_mass_model.setNull();

      } else if (mlower == "constant") {

        Real const_coeff(0.5);
        pp_drag.query("virtual_mass.constant", const_coeff);

        m_virtual_mass_model.setConstant(const_coeff);

      } else if (mlower == "zuber") {

        m_virtual_mass_model.setZuber();

      } else if (mlower == "nijssen") {

        m_virtual_mass_model.setNijssen();

      } else if (mlower != "none") {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Unknown virtual mass model: " << virtual_mass_model_str << "\n"
          << "Please correct the input deck.";
        return 1;
      }
    }
    if (include_virtual_mass()) {
      reporter::Log(reporter::Status) << "Virtual mass model: " << virtual_mass_model_str;
    }
  }

  // Convection model type
  if (a_couple_energy) {

    ParmParse pp_conv("mfix.convection");

    std::string convection_model_str = "RanzMarshall";
    if(!pp_conv.queryAdd("model", convection_model_str)) {
      reporter::Log(reporter::Status)
        << "Convection model not specified. Setting default:\n"
        << "  mfix.convection.model = " << convection_model_str;
    }
    std::string mlower(toLower(convection_model_str));

    if (mlower == "ranzmarshall") {

      m_convection_model.setRanzMarshall();

    } else if (mlower == "gunn") {

      m_convection_model.setGunn();

    } else if (mlower == "nullconvection" || mlower == "none") {

      m_convection_model.setNullConvection();

    } else {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Unknown convection model: " << convection_model_str << "\n"
        << "Please correct the input deck.";
      return 1;
    }
  }

  ParmParse pp_depop("mfix.deposition");

  { std::string filter_str = "none";
    pp_depop.query("filter", filter_str);
    std::string flower(toLower(filter_str));

    Real filter_size = -1;

    if (flower == "constant") {

      if (pp_depop.query("filter.constant", filter_size) ) {

        filter_size = (filter_size*filter_size)/(16.0*std::log(2.0));
        m_depop_filter.setConstantSize(filter_size);

      } else {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Constant size deposition filter requires a size.\n"
          << "Specify the following entry in the inputs file:\n"
          << "   mfix.deposition.filter.constant = amrex::Real";
        return 1;
      }

    } else if (flower == "variable") {

      if (pp_depop.query("filter.variable", filter_size) ) {

        m_depop_filter.setVariableSize(filter_size);

        amrex::Real mineps;
        pp_depop.query("filter.variable.min_eps", mineps);

        if (mineps > 0. && mineps <= 1.) {

          m_depop_filter.setMinSolidsVolfrac(mineps);

        } else {

          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Variable deposition filter - invalid minimum solids volume fraction.\n"
            << "The value must be in the interval (0,1]. Correct the inputs file:\n"
            << "   mfix.deposition.filter.variable.min_eps";
          return 1;
        }

      } else {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Variable deposition filter requires a sample size.\n"
          << "Specify the following entry in the inputs file:\n"
          << "   mfix.deposition.filter.variable = amrex::Real";
        return 1;
      }
    }

    if ( FilterDeposition() ) {
      reporter::Log(reporter::Status)
        << "Deposition filtering is enabled";
    }
  }

  return 0;
}

int CouplingOp::
reset ()
{
  m_on_same_grids.clear();

  for (auto& lev_interp  : m_interp )
  { if (lev_interp  != nullptr) { delete lev_interp;  } }
  m_interp.clear();

  return 0;
}

void CouplingOp::
setup ( Vector<Geometry>& a_geom,
        Vector< MultiFab* > const& a_epf,
        Vector< MultiFab* > const& a_rho,
        Vector< MultiFab* > const& a_vel,
        Vector< MultiFab* > const& a_Tf,
        Vector< MultiFab* > const& a_hf,
        Vector< MultiFab* > const& a_Xfk,
        Vector<EBFArrayBoxFactory const*> const& a_ebfactory)
{
  int const nlev( a_ebfactory.size() );

  m_on_same_grids.resize(nlev,1);

  m_interp.resize(nlev);

  for ( int lev(0); lev<nlev; ++lev) {

    const DistributionMapping& fluid_dm = a_epf[lev]->DistributionMap();
    const BoxArray&            fluid_ba = a_epf[lev]->boxArray();

    const DistributionMapping& eb_dm = a_ebfactory[lev]->DistributionMap();
    const BoxArray&            eb_ba = a_ebfactory[lev]->boxArray();

    m_on_same_grids[lev] = ((fluid_dm == eb_dm) && (fluid_ba == eb_ba));

    if (m_on_same_grids[lev]) {
      m_interp[lev] = new MultiFab(fluid_ba, fluid_dm, interp_comps(),
          interp_ng(), MFInfo(), a_rho[lev]->Factory());
    } else {
      m_interp[lev] = new MultiFab(eb_ba, eb_dm, interp_comps(),
          interp_ng(), MFInfo(), *a_ebfactory[lev]);
    }

#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
    Real const covered_val = 1.e40;

    EB_set_covered(*a_vel[lev], 0, 3, 1, covered_val);
    EB_set_covered(*a_epf[lev], 0, 1, 1, covered_val);
    EB_set_covered(*a_rho[lev], 0, 1, 1, covered_val);

    if ( include_energy() ) {
      EB_set_covered(*a_Tf[lev], 0, 1, 1, covered_val);
      EB_set_covered(*a_hf[lev], 0, 1, 1, covered_val);
    }

    if ( include_species() ) {
      EB_set_covered(*a_Xfk[lev], 0, nspecies(), 1, covered_val);
    }
#endif

    if (OnSameGrids(lev)) {

      // Copy fluid velocity
      MultiFab::Copy(*m_interp[lev], *a_vel[lev], 0, m_interp_idxs->vel_g, 3, interp_ng());

      // Copy volume fraction
      MultiFab::Copy(*m_interp[lev], *a_epf[lev], 0, m_interp_idxs->ep_g, 1, interp_ng());

      // Copy fluid density
      MultiFab::Copy(*m_interp[lev], *a_rho[lev], 0, m_interp_idxs->ro_g, 1, interp_ng());

      if ( include_energy() ) {
        // Copy fluid temperature
        MultiFab::Copy(*m_interp[lev], *a_Tf[lev],  0, m_interp_idxs->T_g, 1, interp_ng());
        MultiFab::Copy(*m_interp[lev], *a_hf[lev],  0, m_interp_idxs->h_g, 1, interp_ng());
      }

      if ( include_species() ) {
        // Copy fluid species
        MultiFab::Copy(*m_interp[lev], *a_Xfk[lev], 0, m_interp_idxs->X_gk,
          m_interp_idxs->nspecies, interp_ng());
      }

    } else {

      // Copy fluid velocity
      m_interp[lev]->ParallelCopy(*a_vel[lev], 0, m_interp_idxs->vel_g, 3, interp_ng(), interp_ng());

      // Copy volume fraction
      m_interp[lev]->ParallelCopy(*a_epf[lev], 0, m_interp_idxs->ep_g, 1, interp_ng(), interp_ng());

      // Copy fluid density
      m_interp[lev]->ParallelCopy(*a_rho[lev], 0, m_interp_idxs->ro_g, 1, interp_ng(), interp_ng());

      if ( include_energy() ) {
        // Copy fluid temperature
        m_interp[lev]->ParallelCopy(*a_Tf[lev],  0, m_interp_idxs->T_g, 1, interp_ng(), interp_ng());
        m_interp[lev]->ParallelCopy(*a_hf[lev],  0, m_interp_idxs->h_g, 1, interp_ng(), interp_ng());
      }

      if ( include_species() ) {
        // Copy fluid species
        m_interp[lev]->ParallelCopy(*a_Xfk[lev], 0, m_interp_idxs->X_gk,
          m_interp_idxs->nspecies, interp_ng(), interp_ng());
      }
    }

    // FillBoundary on interpolation MultiFab
    m_interp[lev]->FillBoundary(a_geom[lev].periodicity());

  } // lev

}
