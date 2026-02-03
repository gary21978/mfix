#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <mfix_rw.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

#ifdef AMREX_USE_ASCENT
#include <AMReX_Conduit_Blueprint.H>
#include <ascent.hpp>
#endif

using namespace amrex;

void
MFIXReadWrite::WriteAscentFile (int nstep, const Real time) const
{
#ifdef AMREX_USE_ASCENT
  BL_PROFILE("mfix::WriteAscentFile()");

  //amrex::Print() << "Writing Ascent output\n";

  if (!m_ascent_actions_yaml.empty()) {

    conduit::Node node;

    Vector<std::string> pltFldNames;
    Vector< MultiFab* > mf(nlev);

    Vector<int> level_steps(nlev);

    int const ngrow(0);
    int ncomp = 1;

    pltFldNames.push_back("volfrac");

    if (fluid.solve()) {

      ncomp += 4;

      pltFldNames.push_back("u_g");
      pltFldNames.push_back("v_g");
      pltFldNames.push_back("w_g");

      pltFldNames.push_back("ep_g");

      // Temperature in fluid
      if (fluid.solve_enthalpy()) {
        pltFldNames.push_back("T_g");
        ncomp += 1;
      }

      if ( fluid.solve_species()) {
        for (std::string specie: fluid.species_names()) {
          pltFldNames.push_back("Xg_"+specie);
          ncomp += 1;
        }
      }

      AMREX_ALWAYS_ASSERT(pltFldNames.size() == ncomp);

    }


    for (int lev(0); lev < nlev; ++lev) {

      level_steps[lev] = nstep;

      mf[lev] = new MultiFab(grids[lev], dmap[lev], ncomp, ngrow,  MFInfo(), *ebfactory[lev]);

      MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, 0, 1, 0);

      if (fluid.solve()) {

        MultiFab::Copy(*mf[lev], *(leveldata_const().vel_const(lev)), 0, 1, 1, 0);
        MultiFab::Copy(*mf[lev], *(leveldata_const().vel_const(lev)), 1, 2, 1, 0);
        MultiFab::Copy(*mf[lev], *(leveldata_const().vel_const(lev)), 2, 3, 1, 0);
        MultiFab::Copy(*mf[lev], *(leveldata_const().epf_const(lev)), 0, 4, 1, 0);

        int lc=5;
        if (fluid.solve_enthalpy()) {
          MultiFab::Copy(*mf[lev], *(leveldata_const().T_const(lev)), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid species mass fractions
        if (fluid.solve_species()) {
          MultiFab::Copy(*mf[lev], *(leveldata_const().X_const(lev)), 0, lc, fluid.nspecies(), 0);
          lc += fluid.nspecies();
        }
      }

      amrex::EB_set_covered(*mf[lev], 0.0);

      amrex::MultiLevelToBlueprint(nlev, amrex::GetVecOfConstPtrs(mf),
          pltFldNames, geom, time, level_steps, ref_ratio, node);

    }


    if ( m_dem.solve() || m_pic.solve() ) {

      Vector<std::string> real_comp_names;
      Vector<std::string>  int_comp_names;

      real_comp_names.push_back("radius");
      real_comp_names.push_back("density");

      real_comp_names.push_back("velx");
      real_comp_names.push_back("vely");
      real_comp_names.push_back("velz");

      if(m_dem.solve()){
        real_comp_names.push_back("omegax");
        real_comp_names.push_back("omegay");
        real_comp_names.push_back("omegaz");
      } else {
        real_comp_names.push_back("grad_tau_x");
        real_comp_names.push_back("grad_tau_y");
        real_comp_names.push_back("grad_tau_z");
      }

      real_comp_names.push_back("statwt");
      real_comp_names.push_back("dragcoeff");
      real_comp_names.push_back("dragx");
      real_comp_names.push_back("dragy");
      real_comp_names.push_back("dragz");

      real_comp_names.push_back("temperature");

      // Currently runtime particle variables are not supported
      // by AMReX.
      int_comp_names.push_back("phase");
      int_comp_names.push_back("state");
#if MFIX_POLYDISPERSE
      int_comp_names.push_back("ptype");
#endif

      amrex::ParticleContainerToBlueprint(*pc,
                real_comp_names, int_comp_names, node);
    }



    // for the MPI case, provide the mpi comm
    ascent::Ascent ascent;

    conduit::Node opts;
    opts["exceptions"] = "catch";
    opts["actions_file"] = m_ascent_actions_yaml;
    opts["mpi_comm"] = MPI_Comm_c2f(ParallelDescriptor::Communicator());

    ascent.open(opts);
    ascent.publish(node);
    conduit::Node actions;
    ascent.execute(actions);
    ParallelDescriptor::Barrier();
    ascent.close();

    for (int lev(0); lev < nlev; ++lev) {
      if(mf[lev] != nullptr) delete mf[lev];
    }
  }

#else
  amrex::ignore_unused(nstep);
  amrex::ignore_unused(time);
#endif
}
