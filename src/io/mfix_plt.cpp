#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_ParmParse.H>
#include <AMReX_EBFArrayBox.H>

#include <mfix.H>
#include <mfix_rw.H>
#include <mfix_pc.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

using namespace amrex;

void
MFIXReadWrite::InitIOPltData ()
{
  if (ooo_debug) amrex::Print() << "InitIOPltData" << std::endl;

  // Variables to simplify checkpoint IO

  pltVarCount = 0;

  ParmParse pp("mfix");

  int plt_ccse_regtest = 0;
  pp.query("plt_regtest", plt_ccse_regtest);

  if (fluid.solve()) {

      pp.query("plt_vel_g",     plt_vel      );
      pp.query("plt_ep_g",      plt_epf      );
      pp.query("plt_p_g",       plt_p_g      );
      pp.query("plt_ro_g",      plt_ro_g     );
      pp.query("plt_MW_g",      plt_MW_g     );
      pp.query("plt_h_g",       plt_h_g      );
      pp.query("plt_T_g",       plt_T_g      );
      pp.query("plt_trac",      plt_trac     );
      pp.query("plt_cp_g",      plt_cp_g     );
      pp.query("plt_k_g",       plt_k_g      );
      pp.query("plt_mu_g",      plt_mu_g     );
      pp.query("plt_vort",      plt_vort     );
      pp.query("plt_volfrac",   plt_volfrac  );
      pp.query("plt_gradp_g",   plt_gradp_g  );
      pp.query("plt_X_g",       plt_X_gk     );
      pp.query("plt_D_g",       plt_D_gk     );
      pp.query("plt_cp_gk",     plt_cp_gk    );
      pp.query("plt_h_gk",      plt_h_gk     );
      pp.query("plt_txfr",      plt_txfr     );
      pp.query("plt_chem_txfr", plt_chem_txfr);
      pp.query("plt_proc",      plt_proc     );
      pp.query("plt_proc_p",    plt_proc_p   );

      // Special test for CCSE regression test. Override all individual
      // flags and save all data to plot file.

      if (plt_ccse_regtest != 0) {
        plt_vel       = 1;
        plt_epf       = 1;
        plt_p_g       = 0;
        plt_ro_g      = 1;
        plt_MW_g      = reactions.solve();
        plt_h_g       = 1;
        plt_T_g       = 1;
        plt_trac      = fluid.solve_tracer();
        plt_cp_g      = 1;
        plt_k_g       = 1;
        plt_mu_g      = 1;
        plt_vort      = 1;
        plt_volfrac   = 1;
        plt_gradp_g   = 1;
        plt_X_gk      = fluid.solve_species();
        plt_D_gk      = fluid.solve_species();
        plt_cp_gk     = 0; //fluid.solve_species() )&&)&& fluid.solve_enthalpy();
        plt_h_gk      = 0; //fluid.solve_species() )&&)&& fluid.solve_enthalpy();
        plt_txfr      = 0;
        plt_chem_txfr = fluid.solve_species() && reactions.solve();
        plt_proc      = 0;
        plt_proc_p    = 0;
      }

      // Count the number of variables to save.
      if (plt_vel      == 1) pltVarCount += 3;
      if (plt_gradp_g  == 1) pltVarCount += 3;
      if (plt_epf      == 1) pltVarCount += 1;
      if (plt_p_g      == 1) pltVarCount += 1;
      if (plt_ro_g     == 1) pltVarCount += 1;
      if (plt_MW_g     == 1) pltVarCount += 1;
      if (plt_trac     == 1) pltVarCount += 1;
      if (plt_mu_g     == 1) pltVarCount += 1;
      if (plt_vort     == 1) pltVarCount += 1;
      if (plt_volfrac  == 1) pltVarCount += 1;
      if (plt_proc     == 1) pltVarCount += 1;
      if (plt_proc_p   == 1) pltVarCount += 1;

      if (fluid.solve_enthalpy()) {
        if (plt_T_g  == 1) pltVarCount += 1;
        if (plt_cp_g == 1) pltVarCount += 1;
        if (plt_k_g  == 1) pltVarCount += 1;
        if (plt_h_g  == 1) pltVarCount += 1;
      }

      if (fluid.solve_species()) {
        if (plt_X_gk == 1)  pltVarCount += fluid.nspecies();
        if (plt_D_gk == 1)  pltVarCount += fluid.nspecies();

        if (fluid.solve_enthalpy()) {
          if (plt_cp_gk == 1) pltVarCount += fluid.nspecies();
          if (plt_h_gk == 1)  pltVarCount += fluid.nspecies();
        }
      }

      InterphaseTxfrIndexes txfr_idxs;

      if (m_dem.solve() || m_pic.solve()) {
        if (plt_txfr == 1) pltVarCount += txfr_idxs.count;
      }

      InterphaseChemTxfrIndexes chem_txfr_idxs(fluid.nspecies(), reactions.solve());

      if (fluid.solve_species() && reactions.solve()) {
        if (plt_chem_txfr == 1) pltVarCount += chem_txfr_idxs.count;
      }
    }

}


void
MFIXReadWrite::
GetSolidsIOPltFlags (std::string const solids_region_prefix,
                     Vector<int>& a_write_real_comp_out,
                     Vector<int>& a_write_int_comp_out)
{
  ParmParse pp(solids_region_prefix);

  int plt_ccse_regtest = 0;
  pp.query("plt_regtest", plt_ccse_regtest);

  runtimeRealData const* const rtData_ptr = &(pc->m_runtimeRealData);
  runtimeIntData const* const itData_ptr = &(pc->m_runtimeIntData);
  const int intsize = SoAintData::count + itData_ptr->count;
  a_write_int_comp_out.resize(intsize, 1);

  // Runtime-added variables
  const int size = SoArealData::count + rtData_ptr->count;
  a_write_real_comp_out.resize(size, 1);

  // All flags are true by default so we only need to turn off the
  // variables we don't want if not doing CCSE regression tests.
  if (plt_ccse_regtest == 0) {

    int input_value = 0;
    pp.query("plt_radius", input_value);
    a_write_real_comp_out[SoArealData::radius] = input_value;

    input_value = 0;
    pp.query("plt_ro_p", input_value);
    a_write_real_comp_out[SoArealData::density] = input_value;

    input_value = 1;
    pp.query("plt_vel_p", input_value);
    a_write_real_comp_out[SoArealData::velx] = input_value;
    a_write_real_comp_out[SoArealData::vely] = input_value;
    a_write_real_comp_out[SoArealData::velz] = input_value;

    input_value = 0;
    pp.query("plt_omega_p", input_value);
    a_write_real_comp_out[SoArealData::omegax] = input_value;
    a_write_real_comp_out[SoArealData::omegay] = input_value;
    a_write_real_comp_out[SoArealData::omegaz] = input_value;

    input_value = 0;
    pp.query("plt_statwt", input_value);
    a_write_real_comp_out[SoArealData::statwt] = input_value;

    input_value = 0;
    pp.query("plt_drag_p", input_value);
    a_write_real_comp_out[SoArealData::drag_coeff] = input_value;
    a_write_real_comp_out[SoArealData::vel_source_x] = input_value;
    a_write_real_comp_out[SoArealData::vel_source_y] = input_value;
    a_write_real_comp_out[SoArealData::vel_source_z] = input_value;

    input_value = 0;
    pp.query("plt_T_p", input_value);
    a_write_real_comp_out[SoArealData::temperature] = input_value;  // temperature

    // Int data
    input_value = 0;
    pp.query("plt_phase", input_value);
    a_write_int_comp_out[SoAintData::phase] = input_value;

    input_value = 0;
    pp.query("plt_state", input_value);
    a_write_int_comp_out[SoAintData::state] = input_value;

  } else {

    a_write_real_comp_out[SoArealData::density] = 0;
    a_write_real_comp_out[SoArealData::statwt] = 0;
  }

  if (rtData_ptr->count > 0) {

    int input_value = 0;

    Vector<int> tmp_write_comp(rtData_ptr->count, 0);

    if (rtData_ptr->contains_species()) {
      AMREX_ALWAYS_ASSERT(solids.solve_species());
      input_value = plt_ccse_regtest;
      pp.query("plt_X_s", input_value);
      for(int n(0); n < solids.nspecies(); ++n) {
        tmp_write_comp[rtData_ptr->X_sn+n] = input_value;
      }
    }

    if (rtData_ptr->contains_eps()) {
      AMREX_ALWAYS_ASSERT(m_pic.solve());
      input_value = 0;
      pp.query("plt_ep_s", input_value);
      tmp_write_comp[rtData_ptr->ep_s] = input_value;
    }

    if (rtData_ptr->contains_vm()) {
      input_value = 0;
      pp.query("plt_vm_coeff", input_value);
      tmp_write_comp[rtData_ptr->vm_coeff] = input_value;
    }

    if (rtData_ptr->contains_acc()) {
      input_value = 0;
      pp.query("plt_acceleration", input_value);
      tmp_write_comp[rtData_ptr->acceleration  ] = input_value;
      tmp_write_comp[rtData_ptr->acceleration+1] = input_value;
      tmp_write_comp[rtData_ptr->acceleration+2] = input_value;
    }

    int const count_tan_history = 3*(rtData_ptr->max_contacts);
    if (count_tan_history > 0) {
       input_value = 0;
       pp.query("plt_pft_neighbor", input_value);
       for(int n(0); n < count_tan_history; ++n) {
          tmp_write_comp[rtData_ptr->pft_neighbor_idx+n] = input_value;
       }
    }

    for (int idx(0); idx<rtData_ptr->count; ++idx) {
      a_write_real_comp_out[SoArealData::count + idx] = tmp_write_comp[idx];
    }
  }


#if MFIX_POLYDISPERSE
  int input_value0 = 0;
  pp.query("plt_ptype", input_value0);
  a_write_int_comp_out[SoAintData::ptype] = input_value0;
#endif

  if (itData_ptr->count > 0) {
    int input_value = 0;
    pp.query("plt_pft_neighbor_flags", input_value);
    const int start = SoAintData::count + itData_ptr->cpu_id_idx;
    for(int n(0); n < itData_ptr->count; ++n) {
      a_write_int_comp_out[start+n] = input_value;
    }
  }

}

void
MFIXReadWrite::
GetSolidsNames () {

  // Vectors of names for solids plot
  real_comp_names.clear();

  real_comp_names.push_back("radius");
  real_comp_names.push_back("density");

  real_comp_names.push_back("velx");
  real_comp_names.push_back("vely");
  real_comp_names.push_back("velz");

  if (m_dem.solve()){
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

  runtimeRealData const* const rtData_ptr = &(pc->m_runtimeRealData);

  if (rtData_ptr->count > 0) {

    Vector<std::string> tmp_comp_name(rtData_ptr->count);

    if (rtData_ptr->contains_species()) {
      AMREX_ALWAYS_ASSERT(solids.solve_species());
      for(int n(0); n < solids.nspecies(); ++n) {
        tmp_comp_name[rtData_ptr->X_sn+n] = "X_"+solids.species_names(n);
      }
    }

    if (rtData_ptr->contains_eps()) {
      AMREX_ALWAYS_ASSERT(m_pic.solve());
      tmp_comp_name[rtData_ptr->ep_s] = "ep_s";
    }

    if (rtData_ptr->contains_vm()) {
      tmp_comp_name[rtData_ptr->vm_coeff] = "vm_coeff";
    }

    if (rtData_ptr->contains_acc()) {
      tmp_comp_name[rtData_ptr->acceleration  ] = "acc_x";
      tmp_comp_name[rtData_ptr->acceleration+1] = "acc_y";
      tmp_comp_name[rtData_ptr->acceleration+2] = "acc_z";
    }

    int const count_tan_history = 3*(rtData_ptr->max_contacts);

    if (count_tan_history > 0) {
       for(int n(0); n < count_tan_history; ++n) {
          tmp_comp_name[rtData_ptr->pft_neighbor_idx+n] =
            "PFT_NEIGHBOR_"+std::to_string(n);
       }
    }

    for (int idx(0); idx<rtData_ptr->count; ++idx) {
      real_comp_names.push_back(tmp_comp_name[idx]);
    }
  } // end of runtime real names

  //--------------- Provide names for integer data ----------------//
  int_comp_names.clear();

  int_comp_names.push_back("phase");
  int_comp_names.push_back("state");
#if MFIX_POLYDISPERSE
  int_comp_names.push_back("ptype");
#endif

  runtimeIntData const* const itData_ptr = &(pc->m_runtimeIntData);

  // Tangential history
  for (int n_th(0); n_th < itData_ptr->count; ++n_th) {
    int_comp_names.push_back("PFT_NEIGHBOR_FLAG_"+std::to_string(n_th));
  }
}




void
MFIXReadWrite::WritePlotFile (std::string& plot_file_in, int nstep, Real time)
{
    // If we've already written this plotfile, don't do it again!
    if (nstep == last_plt) return;

    // Now set last_plt to nstep ...
    last_plt = nstep;

    BL_PROFILE("mfix::WritePlotFile()");

    const std::string& plotfilename = amrex::Concatenate(
        plot_file_in, nstep, m_min_digits
        );

    amrex::Print() << "  Writing plotfile " << plotfilename <<  " at time " << time << std::endl;

    if (pltVarCount > 0) {

      const int ngrow = 0;

      Vector<std::string> pltFldNames;
      Vector< std::unique_ptr<MultiFab> > mf(nlev);

      // Velocity components
      if (plt_vel     == 1) {
        pltFldNames.push_back("u_g");
        pltFldNames.push_back("v_g");
        pltFldNames.push_back("w_g");
      }

      // Pressure gradient
      if (plt_gradp_g == 1) {
        pltFldNames.push_back("gpx");
        pltFldNames.push_back("gpy");
        pltFldNames.push_back("gpz");
      }

      // Fluid volume fraction
      if (plt_epf == 1)
        pltFldNames.push_back("ep_g");

      // Fluid pressure
      if (plt_p_g == 1)
        pltFldNames.push_back("p_g");

      // Fluid density
      if (plt_ro_g == 1)
        pltFldNames.push_back("ro_g");

      // Fluid molecular weight
      if (plt_MW_g == 1)
        pltFldNames.push_back("MW_g");

      // Fluid enthalpy
      if (fluid.solve_enthalpy() && plt_h_g == 1)
        pltFldNames.push_back("h_g");

      // Temperature in fluid
      if (fluid.solve_enthalpy() && plt_T_g == 1)
        pltFldNames.push_back("T_g");

      // Tracer in fluid
      if (plt_trac == 1)
        pltFldNames.push_back("trac");

      // Specific heat
      if (fluid.solve_enthalpy() && plt_cp_g == 1)
        pltFldNames.push_back("cp_g");

      // Thermal conductivity
      if (fluid.solve_enthalpy() && plt_k_g == 1)
        pltFldNames.push_back("k_g");

      // Fluid molecular viscosity
      if (plt_mu_g == 1)
        pltFldNames.push_back("mu_g");

      // vorticity
      if (plt_vort == 1)
        pltFldNames.push_back("vort");

      // EB cell volume fraction
      if (plt_volfrac == 1)
        pltFldNames.push_back("volfrac");

      // rank of fluid grids
      if (plt_proc == 1)
        pltFldNames.push_back("proc");

      // rank of particle grids
      if (plt_proc_p == 1)
        pltFldNames.push_back("proc_p");

      // Fluid species mass fractions
      if (fluid.solve_species() && plt_X_gk == 1)
        for (std::string specie: fluid.species_names())
          pltFldNames.push_back("X_"+specie+"_g");

      // Fluid species mass diffusivities
      if (fluid.solve_species() && plt_D_gk == 1)
        for (std::string specie: fluid.species_names())
          pltFldNames.push_back("D_"+specie+"_g");

      // Fluid species specific heat
      if (fluid.solve_species() && fluid.solve_enthalpy() && plt_cp_gk == 1)
        for (std::string specie: fluid.species_names())
          pltFldNames.push_back("cp_"+specie+"_g");

      // Fluid species enthalpy
      if (fluid.solve_species() && fluid.solve_enthalpy() && plt_h_gk == 1)
        for (std::string specie: fluid.species_names())
          pltFldNames.push_back("h_"+specie+"_g");

      // Fluid species density reaction rates
      if (plt_txfr == 1) {
        pltFldNames.push_back("drag_x");
        pltFldNames.push_back("drag_y");
        pltFldNames.push_back("drag_z");
        pltFldNames.push_back("beta");
        pltFldNames.push_back("gammaTp");
        pltFldNames.push_back("gamma");
      }

      // Fluid species density reaction rates
      if (fluid.solve_species() && reactions.solve() && plt_chem_txfr == 1) {
        for(std::string specie: fluid.species_names())
          pltFldNames.push_back("chem_ro_txfr_"+specie);

        pltFldNames.push_back("chem_h_txfr");
      }

      for (int lev = 0;  lev < nlev; ++lev)
      {
        // Multifab to hold all the variables -- there can be only one!!!!
        const int ncomp = pltVarCount;
        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], ncomp, ngrow,  MFInfo(), *ebfactory[lev]);

        int lc=0;

        // Velocity components
        if (plt_vel   == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().vel(lev)), 0, lc, 3, 0);
          lc += 3;
        }

        // Pressure gradient
        if (plt_gradp_g == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().grad_p(lev)), 0, lc, 3, 0);
          lc += 3;
        }

        // Fluid volume fraction
        if (plt_epf == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().epf(lev)), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid pressure
        if (plt_p_g == 1) {
          amrex::average_node_to_cellcenter(*mf[lev], lc, *(leveldata().pert_p(lev)), 0, 1);
          lc += 1;
        }

        // Fluid density
        if (plt_ro_g == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().rho(lev)), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid molecular weight
        if (plt_MW_g == 1) {

          const int nspecies_g = fluid.nspecies();

          const auto& fluid_parms = fluid.parameters<run_on>();
          const int fluid_is_a_mixture = fluid.isMixture();

          MultiFab& epf = *(leveldata().epf(lev));

          MultiFab MW_g(epf.boxArray(), epf.DistributionMap(), epf.nComp(),
                        epf.nGrow(), MFInfo(), epf.Factory());

          for (MFIter mfi(epf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& MW_g_array = MW_g.array(mfi);
            Array4<Real const> const& X_gk_array = leveldata().X_const(lev,mfi);

            ParallelFor(bx, [MW_g_array,X_gk_array,nspecies_g,fluid_is_a_mixture,
                fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              if (fluid_is_a_mixture) {
                Real MW_g_loc(0);

                for (int n(0); n < nspecies_g; ++n) {
                  const Real MW_gk = fluid_parms.get_MW_gk(n);

                  MW_g_loc += X_gk_array(i,j,k,n) / MW_gk;
                }

                MW_g_array(i,j,k) = 1. / MW_g_loc;
              } else {
                MW_g_array(i,j,k) = fluid_parms.get_MW_g();
              }
            });
          }

          EB_set_covered(MW_g, 0, MW_g.nComp(), MW_g.nGrow(), mfix::covered_val);

          MultiFab::Copy(*mf[lev], MW_g, 0, lc, 1, 0);

          lc += 1;
        }

        // Fluid enthalpy
        if (fluid.solve_enthalpy() && plt_h_g == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().h(lev)), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid temperature
        if (fluid.solve_enthalpy() && plt_T_g == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().T(lev)), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid tracer
        if (plt_trac == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().tracer(lev)), 0, lc, 1, 0);
          lc += 1;
        }

        // Specific heat
        if (fluid.solve_enthalpy() && plt_cp_g == 1) {

          const auto fluid_props = fluid.props.data<run_on>();

          MultiFab& Tf = *(leveldata().T(lev));

          MultiFab cp_g(Tf.boxArray(), Tf.DistributionMap(), Tf.nComp(),
                        Tf.nGrow(), MFInfo(), Tf.Factory());

          for (MFIter mfi(Tf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& cp_g_array = cp_g.array(mfi);
            Array4<Real const> const& Tf_array  = Tf.const_array(mfi);

            Array4<Real const> const& X_gk_array = leveldata().X_const(lev,mfi);

            ParallelFor(bx, [cp_g_array,Tf_array,X_gk_array,fluid_props]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              cp_g_array(i,j,k) = fluid_props.specificHeat(IntVect(i,j,k), Tf_array, X_gk_array);
            });
          }

          EB_set_covered(cp_g, 0, cp_g.nComp(), cp_g.nGrow(), mfix::covered_val);

          MultiFab::Copy(*mf[lev], cp_g, 0, lc, 1, 0);
          lc += 1;
        }

        // Thermal conductivity
        if (fluid.solve_enthalpy() && plt_k_g == 1) {

          const auto& fluid_parms = fluid.parameters<run_on>();

          MultiFab& Tf = *(leveldata().T(lev));

          MultiFab k_g(Tf.boxArray(), Tf.DistributionMap(), Tf.nComp(),
                       Tf.nGrow(), MFInfo(), Tf.Factory());

          for (MFIter mfi(Tf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& k_g_array = k_g.array(mfi);
            Array4<Real const> const& Tf_array = Tf.const_array(mfi);

            ParallelFor(bx, [k_g_array,Tf_array,fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              k_g_array(i,j,k) = fluid_parms.calc_k_g(Tf_array(i,j,k));
            });
          }

          EB_set_covered(k_g, 0, k_g.nComp(), k_g.nGrow(), mfix::covered_val);

          MultiFab::Copy(*mf[lev], k_g, 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid molecular viscosity
        if (plt_mu_g == 1) {

          MultiFab& epf = *(leveldata().epf(lev));

          MultiFab mu_g(epf.boxArray(), epf.DistributionMap(), epf.nComp(),
                        epf.nGrow(), MFInfo(), epf.Factory());

          for (MFIter mfi(epf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& mu_g_array = mu_g.array(mfi);
            Array4<Real const> const& Tf_array   = leveldata().T_const(lev,mfi);
            Array4<Real const> const& X_gk_array = leveldata().X_const(lev,mfi);

            const auto fluid_props = fluid.props.data<run_on>();

            ParallelFor(bx, [mu_g_array,Tf_array,X_gk_array,fluid_props]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               mu_g_array(i,j,k) = fluid_props.molViscosity(i, j, k, Tf_array, X_gk_array);
            });
          }

          EB_set_covered(mu_g, 0, mu_g.nComp(), mu_g.nGrow(), mfix::covered_val);

          MultiFab::Copy(*mf[lev], mu_g, 0, lc, 1, 0);
          lc += 1;
        }

        // vorticity
        if (plt_vort == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().vorticity(lev)), 0, lc, 1, 0);
          lc += 1;
        }

        // EB cell volume fraction
        if (plt_volfrac == 1) {
          if (ebfactory[lev]) {
            MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, lc, 1, 0);
          } else {
            mf[lev]->setVal(1.0, lc, 1, 0);
          }
          lc += 1;
        }

        // rank of fluid grids
        if( plt_proc == 1 ) {

          Real const proc( static_cast<Real>(ParallelDescriptor::MyProc()) );
          mf[lev]->setVal(proc, lc, 1, 0);

          lc += 1;
        }

        // rank of particle grids
        if ( plt_proc_p == 1 ) {

          Real const proc( static_cast<Real>(ParallelDescriptor::MyProc()) );

          const DistributionMapping& pc_dm = pc->ParticleDistributionMap(lev);
          const BoxArray&            pc_ba = pc->ParticleBoxArray(lev);

          if ( (grids[lev] == pc_ba) && (dmap[lev] == pc_dm) ) {
            mf[lev]->setVal(proc, lc, 1, 0);

          } else {

            MultiFab pc_procs(pc_ba, pc_dm, 1, 0);
            pc_procs.setVal(proc);
            mf[lev]->ParallelCopy(pc_procs, 0, lc, 1, 0, 0);

          }

          lc += 1;
        }

        // Fluid species mass fractions
        if (fluid.solve_species() && plt_X_gk == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().X(lev)), 0, lc, fluid.nspecies(), 0);
          lc += fluid.nspecies();
        }

        // Species mass fraction
        if (fluid.solve_species() && plt_D_gk == 1) {

          const auto& fluid_parms = fluid.parameters<run_on>();
          const int nspecies_g = fluid.nspecies();

          MultiFab& X_gk = *(leveldata().X(lev));

          MultiFab D_gk(X_gk.boxArray(), X_gk.DistributionMap(), nspecies_g,
                        X_gk.nGrow(), MFInfo(), X_gk.Factory());

          for (MFIter mfi(X_gk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& D_gk_array = D_gk.array(mfi);

            ParallelFor(bx, [D_gk_array,nspecies_g,fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              for (int n(0); n < nspecies_g; ++n) {
                D_gk_array(i,j,k,n) = fluid_parms.get_D_g();
              }
            });
          }

          EB_set_covered(D_gk, 0, D_gk.nComp(), D_gk.nGrow(), mfix::covered_val);

          MultiFab::Copy(*mf[lev], D_gk, 0, lc, D_gk.nComp(), 0);

          lc += D_gk.nComp();
        }

        // Fluid species specific heat
        if (fluid.solve_species() && fluid.solve_enthalpy() && plt_cp_gk == 1) {

          const auto fluid_props = fluid.props.data<run_on>();
          const int nspecies_g = fluid.nspecies();

          MultiFab& X_gk = *(leveldata().X(lev));

          MultiFab cp_gk(X_gk.boxArray(), X_gk.DistributionMap(), X_gk.nComp(),
                         X_gk.nGrow(), MFInfo(), X_gk.Factory());

          for (MFIter mfi(X_gk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& cp_gk_array = cp_gk.array(mfi);
            Array4<Real const> const& Tf_array  = leveldata().T_const(lev,mfi);

            ParallelFor(bx, [cp_gk_array,Tf_array,nspecies_g,fluid_props]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              for (int n(0); n < nspecies_g; ++n) {
                const Real Tf = Tf_array(i,j,k);

                cp_gk_array(i,j,k,n) = fluid_props.specificHeat(n,Tf);
              }
            });
          }

          EB_set_covered(cp_gk, 0, cp_gk.nComp(), cp_gk.nGrow(), mfix::covered_val);

          MultiFab::Copy(*mf[lev], cp_gk, 0, lc, nspecies_g, 0);

          lc += nspecies_g;
        }

        // Fluid species enthalpy
        if (fluid.solve_species() && plt_h_gk == 1) {

          const auto fluid_props = fluid.props.data<run_on>();
          const int nspecies_g = fluid.nspecies();

          MultiFab& X_gk = *(leveldata().X(lev));

          MultiFab h_gk(X_gk.boxArray(), X_gk.DistributionMap(), X_gk.nComp(),
                        X_gk.nGrow(), MFInfo(), X_gk.Factory());

          for (MFIter mfi(X_gk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            const EBFArrayBox& Xgk_fab = static_cast<EBFArrayBox const&>(X_gk[mfi]);
            const EBCellFlagFab& flags = Xgk_fab.getEBCellFlagFab();

            Box const& bx = mfi.tilebox();

            Array4<      Real> const& h_gk_array = h_gk.array(mfi);
            Array4<const Real> const& Tf_array  = leveldata().T_const(lev,mfi);

            auto const& flags_arr = flags.const_array();

            ParallelFor(bx, [h_gk_array,Tf_array,nspecies_g,fluid_props,flags_arr]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

              const Real Tg_loc = Tf_array(i,j,k);

              for (int n(0); n < nspecies_g; ++n)
                h_gk_array(i,j,k,n) = fluid_props.enthalpy(n, Tg_loc, cell_is_covered);
            });
          }

          EB_set_covered(h_gk, 0, h_gk.nComp(), h_gk.nGrow(), mfix::covered_val);

          MultiFab::Copy(*mf[lev], h_gk, 0, lc, nspecies_g, 0);

          lc += nspecies_g;
        }

        InterphaseTxfrIndexes txfr_idxs;

        if (plt_txfr == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().txfr(lev)), 0, lc, txfr_idxs.count, 0);
          lc += txfr_idxs.count;
        }

        InterphaseChemTxfrIndexes chem_txfr_idxs(fluid.nspecies(), reactions.solve());

        if (fluid.solve_species() && reactions.solve() && plt_chem_txfr == 1) {
          MultiFab::Copy(*mf[lev], *(leveldata().chem_txfr(lev)),
            0, lc, chem_txfr_idxs.count, 0);

          lc += chem_txfr_idxs.count;
        }

      }

      // Cleanup places where we have no data.
      Vector<const MultiFab*> mf2(nlev);
      for (int lev = 0;  lev < nlev; ++lev) {
        EB_set_covered(*mf[lev], 0.0);
        mf2[lev] = mf[lev].get();
      }

      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfile(plotfilename, nlev, mf2, pltFldNames,
                                     geom, time, istep, ref_ratio);


      // no fluid
    } else {

      // Some post-processing tools (such as yt) might still need some basic
      // MultiFab header information to function. We provide this here by
      // creating an "empty" plotfile header (which essentially only contains
      // the BoxArray information). Particle data is saved elsewhere.

      Vector< std::unique_ptr<MultiFab> > mf(nlev);
      Vector<std::string>  names;
      // NOTE: leave names vector empty => header should reflect nComp = 0
      //names.insert(names.end(), "placeholder");

      // Create empty MultiFab containing the right BoxArray (NOTE: setting
      // nComp = 1 here to avoid assertion fail in debug build).
      for (int lev = 0; lev < nlev; ++lev)
        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], 1, 0);

      Vector<const MultiFab*> mf2(nlev);

      for (int lev = 0; lev < nlev; ++lev)
        mf2[lev] = mf[lev].get();

      // Write only the Headers corresponding to the "empty" mf/mf2 MultiFabs
      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfileHeaders(plotfilename, nlev, mf2, names,
                                            geom, time, istep, ref_ratio);

    }

    WriteJobInfo(plotfilename);

    if (m_dem.solve() || m_pic.solve()) {

        pc->WritePlotFile(plotfilename, "particles", write_real_comp,
                          write_int_comp, real_comp_names, int_comp_names);
    }
}


void MFIXReadWrite::
WriteSolidsPlotFile ( SolidsPlotRegion& plot_region,
                      int nstep,
                      Real time)
{
  if ((m_dem.solve() || m_pic.solve()) && (solids_plot_regions() == true)) {

    // If we've already written this plotfile, don't do it again!
    if (nstep == plot_region.m_last_solids_plt) return;

    // Now set last_solids_plt to nstep ...
    plot_region.m_last_solids_plt = nstep;

    const RealBox region_extents = plot_region.m_region_extents;
    const std::string& region_name = plot_region.m_region_name;

    std::string basename = amrex::Concatenate("plt",nstep, m_min_digits);
    std::string fname = region_name + "/" + basename;

    amrex::Print() << "  Writing solids plotfile " << fname <<  " at time " << time << '\n';

    const int plot_types_nb = plot_region.m_h_plot_types.size();
    int* plot_types_ptr = plot_region.m_d_plot_types.dataPtr();

    auto F = [region_extents,plot_types_ptr,plot_types_nb]
      AMREX_GPU_DEVICE (const MFIXParticleContainer::SuperParticleType& p,
                        const amrex::RandomEngine&) noexcept -> bool
    {
      int particle_in_region = static_cast<int>(region_extents.contains(p.pos()));

      int type_found = static_cast<int>(plot_types_nb == 0);
      for (int n(0); n < plot_types_nb; ++n) {
        if (p.idata(SoAintData::phase) == plot_types_ptr[n]) {
          type_found = 1;
          break;
        }
      }

      return particle_in_region && type_found;
    };

    amrex::Vector<std::string>& fluid_vars = plot_region.m_plot_fluid_vars;
    const int ncomp = fluid_vars.size();
    const int ngrow = 0;

    Vector<std::unique_ptr<MultiFab>> mf(nlev);

    if (ncomp > 0) {

      for (int lev = 0; lev < nlev; ++lev) {

        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], ncomp, ngrow,  MFInfo(), *ebfactory[lev]);

        int lc(0);

        for (std::string& var: fluid_vars) {

          if (var.compare("ep_g") == 0) {

            MultiFab::Copy(*mf[lev], *(leveldata().epf(lev)), 0, lc, 1, 0);
            lc += 1;

          } else if (var.compare("ro_g") == 0) {

            MultiFab::Copy(*mf[lev], *(leveldata().rho(lev)), 0, lc, 1, 0);
            lc += 1;

          } else if (var.compare("T_g") == 0) {

            MultiFab::Copy(*mf[lev], *(leveldata().T(lev)), 0, lc, 1, 0);
            lc += 1;

          } else if (var.compare("volfrac") == 0) {

            if (ebfactory[lev]) {
              MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, lc, 1, 0);
            } else {
              mf[lev]->setVal(1.0, lc, 1, 0);
            }
            lc += 1;
          } else {
            amrex::Abort("Error: invalid fluid variable in solids region plot");
          }
        }
      }

      // Cleanup places where we have no data.
      Vector<const MultiFab*> mf2(nlev);
      for (int lev = 0;  lev < nlev; ++lev) {
        EB_set_covered(*mf[lev], 0.0);
        mf2[lev] = mf[lev].get();
      }

      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfile(fname, nlev, mf2, fluid_vars,
                                     geom, time, istep, ref_ratio);

    // no fluid
    } else {
      // Some post-processing tools (such as yt) might still need some basic
      // MultiFab header information to function. We provide this here by
      // creating an "empty" plotfile header (which essentially only contains
      // the BoxArray information). Particle data is saved elsewhere.

      Vector<std::string> names;
      // NOTE: leave names vector empty => header should reflect nComp = 0
      //names.insert(names.end(), "placeholder");

      // Create empty MultiFab containing the right BoxArray (NOTE: setting
      // nComp = 1 here to avoid assertion fail in debug build).
      for (int lev = 0; lev < nlev; ++lev)
        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], 1, 0);

      Vector<const MultiFab*> mf2(nlev);

      for (int lev = 0; lev < nlev; ++lev)
        mf2[lev] = mf[lev].get();

      // Write only the Headers corresponding to the "empty" mf/mf2 MultiFabs
      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfileHeaders(fname, nlev, mf2, names,
                                            geom, time, istep, ref_ratio);
    }

    WriteJobInfo(fname);

    pc->WritePlotFile(fname, "particles", plot_region.m_write_real_comp,
                      plot_region.m_write_int_comp, real_comp_names, int_comp_names, F);
  }
}


void
MFIXReadWrite::WriteStaticPlotFileParticleLevelSet (const std::string & plotfilename) const
{
    BL_PROFILE("mfix::WriteStaticPlotFileParticleLevelSet()");

    Print() << "  Writing static quantities " << plotfilename << std::endl;

    /****************************************************************************
     *                                                                          *
     * Static (un-changing variables):                                          *
     *     1. level-set data                                                    *
     *     2. volfrac (from EB) data                                            *
     *                                                                          *
     ***************************************************************************/

    Vector<std::string> static_names = {"level_sets", "volfrac"};
    Vector< Vector< MultiFab const* > > static_vars = { amrex::GetVecOfConstPtrs(level_sets) };

    const int ngrow = 0;
    const int ncomp = static_names.size();


    /****************************************************************************
     *                                                                          *
     * Collect variables together into a single multi-component MultiFab        *
     *                                                                          *
     ***************************************************************************/

    Vector<std::unique_ptr<MultiFab>> mf(nlev);
    Vector<const MultiFab *>          mf_ptr(nlev);

    for (int lev = 0;  lev < nlev; lev++)
    {
        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], ncomp, ngrow,
                                             MFInfo(), *particle_ebfactory[lev]);

        // Don't iterate over all ncomp => last component is for volfrac
        for (int dcomp = 0; dcomp < ncomp - 1; dcomp++)
        {
            const BoxArray nd_ba = amrex::convert(grids[lev], IntVect::TheNodeVector());
            MultiFab mf_loc = MFUtil::regrid(nd_ba, dmap[lev], *static_vars[dcomp][lev], true);
            amrex::average_node_to_cellcenter(* mf[lev], dcomp, mf_loc, 0, 1, ngrow);
        }

        if (ebfactory[lev]) {
            MultiFab::Copy(* mf[lev], ebfactory[lev]->getVolFrac(), 0, ncomp - 1, 1, ngrow);

        } else {
            // setVal (value_type val, int comp, int num_comp, int nghost=0)
            mf[lev]->setVal(1.0, ncomp - 1, 1, ngrow);
        }
    }

    for (int lev = 0;  lev < nlev; ++lev)
    {
        // Don't do this (below) as it zeros out the covered cells...
        // EB_set_covered(* mf[lev], 0.0);
        mf_ptr[lev] = mf[lev].get();
    }

    Real time = 0.;
    Vector<int> istep;
    istep.resize(nlev,0);
    amrex::WriteMultiLevelPlotfile(plotfilename, nlev, mf_ptr, static_names,
                                   geom, time, istep, ref_ratio);

    WriteJobInfo(plotfilename);

    Print() << "  Done writing static quantities " << plotfilename << std::endl;
}

void
MFIXReadWrite::WriteStaticPlotFileEBGeometry (const std::string & plotfilename) const
{
    BL_PROFILE("mfix::WriteStaticPlotFileEBGeometry()");

    Print() << "  Writing static quantities " << plotfilename << std::endl;

    // Static unchanging EB geometry variables
    Vector<std::string> static_names =
      {"volfrac",
       "impf",
       "centroid_x", "centroid_y", "centroid_z",
       "bndryarea",
       "bndrycent_x", "bndrycent_y", "bndrycent_z",
       "bndrynorm_x", "bndrynorm_y", "bndrynorm_z",
       "areafrac_x", "areafrac_y", "areafrac_z",
       "facecent_xy", "facecent_xz", "facecent_yx", "facecent_yz", "facecent_zx", "facecent_zy",
       "edgecent_x", "edgecent_y", "edgecent_z"};

    const int ngrow = 0;
    const int ncomp = static_names.size();


    /****************************************************************************
     *                                                                          *
     * Collect variables together into a single multi-component MultiFab        *
     *                                                                          *
     ***************************************************************************/

    Vector<std::unique_ptr<MultiFab>> mf(nlev);
    Vector<const MultiFab *>          mf_ptr(nlev);

    for (int lev = 0;  lev < nlev; lev++)
    {
        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], ncomp, ngrow,
                                             MFInfo(), *ebfactory[lev]);

        int lc = 0;

        // volfrac
        if (ebfactory[lev]) {
            MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, lc, 1, 0);
        } else {
            mf[lev]->setVal(1.0, lc, 1, 0);
        }
        lc += 1;

        // implicit function or levelset
        if (ebfactory[lev]) {
            MultiFab::Copy(*mf[lev], ebfactory[lev]->getLevelSet(), 0, lc, 1, 0);
        } else {
            mf[lev]->setVal(-1.0, lc, 1, 0);
        }
        lc += 1;

        // centroid
        if (ebfactory[lev]) {
            MultiFab temp_mf = ebfactory[lev]->getCentroid().ToMultiFab(0.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mf, 0, lc, 3, 0);
        } else {
            mf[lev]->setVal(0.0, lc, 3, 0);
        }
        lc += 3;

        // boundary area
        if (ebfactory[lev]) {
            MultiFab temp_mf = ebfactory[lev]->getBndryArea().ToMultiFab(0.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mf, 0, lc, 1, 0);
        } else {
            mf[lev]->setVal(0.0, lc, 1, 0);
        }
        lc += 1;

        // boundary centroid
        if (ebfactory[lev]) {
            MultiFab temp_mf = ebfactory[lev]->getBndryCent().ToMultiFab(-1.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mf, 0, lc, 3, 0);
        } else {
            mf[lev]->setVal(-1.0, lc, 3, 0);
        }
        lc += 3;

        // boundary normal
        if (ebfactory[lev]) {
            MultiFab temp_mf = ebfactory[lev]->getBndryNormal().ToMultiFab(0.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mf, 0, lc, 3, 0);
        } else {
            mf[lev]->setVal(0.0, lc, 3, 0);
        }
        lc += 3;

        // area fraction
        if (ebfactory[lev]) {
            MultiFab temp_mfx = ebfactory[lev]->getAreaFrac()[0]->ToMultiFab(1.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mfx, 0, lc, 1, 0);
            lc += 1;
            MultiFab temp_mfy = ebfactory[lev]->getAreaFrac()[1]->ToMultiFab(1.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mfy, 0, lc, 1, 0);
            lc += 1;
            MultiFab temp_mfz = ebfactory[lev]->getAreaFrac()[2]->ToMultiFab(1.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mfz, 0, lc, 1, 0);
            lc += 1;
        } else {
            mf[lev]->setVal(1.0, lc, 3, 0);
            lc += 3;
        }

        // face centroid
        if (ebfactory[lev]) {
            MultiFab temp_mfx = ebfactory[lev]->getFaceCent()[0]->ToMultiFab(0.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mfx, 0, lc, 2, 0);
            lc += 2;
            MultiFab temp_mfy = ebfactory[lev]->getFaceCent()[1]->ToMultiFab(0.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mfy, 0, lc, 2, 0);
            lc += 2;
            MultiFab temp_mfz = ebfactory[lev]->getFaceCent()[2]->ToMultiFab(0.0, 0.0);
            MultiFab::Copy(*mf[lev], temp_mfz, 0, lc, 2, 0);
            lc += 2;
        } else {
            mf[lev]->setVal(0.0, lc, 6, 0);
            lc += 6;
        }

        // edge centroid
        if (ebfactory[lev]) {
            MultiFab temp_mfx = ebfactory[lev]->getEdgeCent()[0]->ToMultiFab(1.0, -1.0);
            MultiFab::Copy(*mf[lev], temp_mfx, 0, lc, 1, 0);
            lc += 1;
            MultiFab temp_mfy = ebfactory[lev]->getEdgeCent()[1]->ToMultiFab(1.0, -1.0);
            MultiFab::Copy(*mf[lev], temp_mfy, 0, lc, 1, 0);
            lc += 1;
            MultiFab temp_mfz = ebfactory[lev]->getEdgeCent()[2]->ToMultiFab(1.0, -1.0);
            MultiFab::Copy(*mf[lev], temp_mfz, 0, lc, 1, 0);
            lc += 1;
        } else {
            mf[lev]->setVal(1.0, lc, 3, 0);
            lc += 3;
        }
    }

    for (int lev = 0;  lev < nlev; ++lev)
    {
        mf_ptr[lev] = mf[lev].get();
    }

    Real time = 0.;
    Vector<int> istep;
    istep.resize(nlev,0);
    amrex::WriteMultiLevelPlotfile(plotfilename, nlev, mf_ptr, static_names,
                                   geom, time, istep, ref_ratio);

    WriteJobInfo(plotfilename);

    Print() << "  Done writing static quantities " << plotfilename << std::endl;
}
