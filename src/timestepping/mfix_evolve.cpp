#include <AMReX_BCUtil.H>

#include <mfix.H>
#include <mfix_utils.H>

void mfix::Evolve ()
{
  if (m_use_new_depdt_algo) {

    if (advection_type() == AdvectionType::MOL) {

      Evolve_MOL();

    } else {

      Evolve_Godunov();

    }

  } else {

    Evolve_OLD();

  }
}

// This subroutine is the driver for the whole time stepping (fluid + particles )
void mfix::Evolve_OLD ()
{
  BL_PROFILE_REGION_START("mfix::Evolve");

  int const t_fluid(0), t_solids(1), t_coupling(2);
  Vector<Real> timings = {0., 0., 0.};

  ///////////////////////////////////////////////////////////////////////////////
  //                           PREPARE FLUID FOR STEP                          //
  ///////////////////////////////////////////////////////////////////////////////
  MFIXStepData stepData(nlev(), grids, dmap, eb()->factory_const(),
      nghost_state(), nghost_force(), bcs() );

  prepareToStep( stepData );

  ///////////////////////////////////////////////////////////////////////////////
  //                           COMPUTE VOLUME FRACTION                         //
  ///////////////////////////////////////////////////////////////////////////////
  if ( fluid.solve() ) {

    Real const timer_volfrac( ParallelDescriptor::second() );

    if (m_dem.solve() || m_pic.solve()) {

      // Swap old and new epf multifabs
      if ( fluid.constraint.include_depdt() ) {
        for (int lev(0); lev < nlev(); ++lev) {
          std::swap(m_level_data.m_level_data[lev]->m_epf,
                    m_level_data.m_level_data[lev]->m_epf_old);
        }
      }

      deposit_volume_to_grid( m_timer.time() );

      if (m_porous_media.nregions() > 0) {
        m_porous_media.impose_vfrac(leveldata().epf());
      }

      // Set the volume fraction to the correct values in inflow faces.
      for (int lev(0); lev < nlev(); ++lev) {
        bcs().set_epf_bcs(lev, m_timer.time(), leveldata().epf(lev), BCType::ext_dir);
      }

      for (int lev(0); lev < nlev(); ++lev) {

        MultiFab::Copy(*stepData.epf_nph(lev), *(leveldata().epf(lev)),
            0, 0, 1, leveldata().epf(lev)->nGrow());

        // Include depdt in the constraint, therefore compute the 'half-time' volume
        // half-time volume fraction as an average of old and new fields.
        if ( fluid.constraint.include_depdt() ) {

          MultiFab::Add(*stepData.epf_nph(lev), *leveldata().epf_old(lev),
              0, 0, 1, leveldata().epf(lev)->nGrow());

          stepData.epf_nph(lev)->mult(0.5, 0, 1, leveldata().epf(lev)->nGrow());

          if ( !almostEqual(m_timer.dt(), Real(0.)) ) {

              Real const inv_dt( 1_rt/m_timer.dt() );

              MultiFab::LinComb(*stepData.depdt(lev), inv_dt, *leveldata().epf(lev), 0,
                  -inv_dt, *leveldata().epf_old(lev), 0, 0, 1,
                  stepData.depdt(lev)->nGrow());
         }

        // Not including depdt in the constraint, therefore copy the new volume
        // fraction into the old.
        } else {

          MultiFab::Copy(*(leveldata().epf_old(lev)), *(leveldata().epf(lev)),
              0, 0, 1, leveldata().epf(lev)->nGrow());

        }
      }

      //const IntVect min_epg_cell = m_rw->mfix_print_min_epg();
      timings[t_coupling] = (ParallelDescriptor::second() - timer_volfrac);

    }
    else { // no dem or pic

      if (m_porous_media.nregions() > 0) {
        m_porous_media.impose_vfrac(leveldata().epf());
      }
    }
  }


  ///////////////////////////////////////////////////////////////////////////////
  //                 COMPUTE FLUID-SOLIDS COUPLING AND REACTIONS               //
  ///////////////////////////////////////////////////////////////////////////////
  if (fluid.solve() && (m_dem.solve() || m_pic.solve() || reactions.solve()) ) {

    Real const timer_calc_txfr( ParallelDescriptor::second() );

    // We copy the value inside the domain to the outside to avoid unphysical
    // volume fractions.
    for ( int lev(0); lev < nlev(); ++lev) {
      bcs().set_epf_bcs(lev, m_timer.time(), leveldata().epf(lev), BCType::foextrap);
    }

    // Setup the interpolation MultiFab
    m_coupling.setup( Geom(), leveldata().epf(), leveldata().rho(), leveldata().vel(),
        leveldata().T(), leveldata().h(), leveldata().X(), eb()->pc_factory_const());

    // Reset the volume fractions back to the correct values at inflow faces.
    for ( int lev(0); lev < nlev(); ++lev) {
      bcs().set_epf_bcs(lev, m_timer.time(), leveldata().epf(lev), BCType::ext_dir);
    }

    if (m_dem.solve() || m_pic.solve()) {

      m_coupling.calc_transfer_coeffs(fluid, pc, geom[0].data());
      deposit_forces_to_grid(m_timer.dt());
    }

    if (reactions.solve()) {

      m_cmi.update_pointers(&reactions, &solids, &fluid, &m_dem, &m_pic, pc,
          &m_timer, m_rw->report_mass_balance);

      m_cmi.calc_chemistry(geom, leveldata().chem_txfr(), m_coupling.get_interp(),
        therm_p, m_timer.time(), m_timer.dt(), m_verbose);

      m_mass_balance->ComputeMassProduction(m_cmi.fluid_mass_change(),
          m_cmi.solids_mass_change());
    }

    m_coupling.reset();

    timings[t_coupling] += (ParallelDescriptor::second() - timer_calc_txfr);
  }


  ///////////////////////////////////////////////////////////////////////////////
  //                            FLUID PREDICTOR STEP                           //
  ///////////////////////////////////////////////////////////////////////////////
  if (fluid.solve()) {

    Real const timer_predictor( ParallelDescriptor::second() );

    compute_dt(stepData.avg_particle_data_const());

    Print() << "\n   Step " << m_timer.nstep()+1 << ": from old_time " \
            << m_timer.time() << " to new time " << m_timer.new_time()
            << " with dt = " << m_timer.dt() << "\n" << std::endl;

    // NOTE: It is important to do fillpatch on the NEW variables,
    // before swapping pointers.
    for (int lev(0); lev < nlev(); lev++) {

      Real const time = timer().time();
      int const ng = nghost_state();

      bcs().fillpatch(lev, time, BCFillVar::vel, leveldata().vel(), ng);
      bcs().fillpatch(lev, time, BCFillVar::rho, leveldata().rho(), ng);

      if (leveldata(lev)->has_enthalpy()) {
        bcs().fillpatch(lev, time, BCFillVar::h, leveldata().h(), ng);
      }
      if (leveldata(lev)->has_temperature()) {
        bcs().fillpatch(lev, time, BCFillVar::T, leveldata().T(), ng);
      }
      if (leveldata(lev)->has_species()) {
        bcs().fillpatch(lev, time, BCFillVar::X, leveldata().X(), ng);
      }
      if (leveldata(lev)->has_tracer()) {
        bcs().fillpatch(lev, time, BCFillVar::tracer, leveldata().tracer(), ng);
      }
      leveldata(lev)->swapOldAndNew();
    }
    if (m_verbose > 0) { Print() << " Finished fillpatch before predictor\n"; }

    std::swap(therm_p, therm_po);

    apply_predictor(stepData, m_timer.time(), m_timer.dt(),
        m_timer.prev_dt(), /*proj_2_pred=*/true);

    timings[t_fluid] += (ParallelDescriptor::second() - timer_predictor);
  }


  ///////////////////////////////////////////////////////////////////////////////
  //              RE-COMPUTE FLUID-SOLIDS COUPLING SKIP REACTIONS              //
  ///////////////////////////////////////////////////////////////////////////////
  if (advection_type() == AdvectionType::MOL) {

    if (fluid.solve() && (m_dem.solve() || m_pic.solve()) ) {

      Real const timer_calc_txfr( ParallelDescriptor::second() );

      // Copy the value inside the domain to the outside to avoid unphysical
      // volume fractions.
      for (int lev(0); lev < nlev(); ++lev) {
        bcs().set_epf_bcs(lev, m_timer.time(), leveldata().epf(lev), BCType::foextrap);
      }

      // Setup the interpolation MultiFab
      m_coupling.setup( Geom(), leveldata().epf(), leveldata().rho(),
          leveldata().vel(), leveldata().T(), leveldata().h(), leveldata().X(),
          eb()->pc_factory_const());

      // Reset the volume fractions back to the correct values at inflow faces.
      for (int lev(0); lev < nlev(); ++lev) {
        bcs().set_epf_bcs(lev, m_timer.time(), leveldata().epf(lev), BCType::ext_dir);
      }

      m_coupling.calc_transfer_coeffs(fluid, pc, geom[0].data());
      deposit_forces_to_grid(m_timer.dt());

      m_coupling.reset();

      timings[t_coupling] += (ParallelDescriptor::second() - timer_calc_txfr);
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  //                           CORRECTOR STEP                              //
  ///////////////////////////////////////////////////////////////////////////
  if (advection_type() == AdvectionType::MOL) {

    if (fluid.solve()) {

      Real const timer_corrector(ParallelDescriptor::second());

      for (int lev(0); lev < nlev(); lev++) {

        Real const time = timer().time();
        int const ng = nghost_state();

        bcs().fillpatch(lev, time, BCFillVar::vel, leveldata().vel(), ng);
        bcs().fillpatch(lev, time, BCFillVar::rho, leveldata().rho(), ng);

        if (leveldata(lev)->has_enthalpy()) {
          bcs().fillpatch(lev, time, BCFillVar::h, leveldata().h(), ng);
        }
        if (leveldata(lev)->has_temperature()) {
          bcs().fillpatch(lev, time, BCFillVar::T, leveldata().T(), ng);
        }
        if (leveldata(lev)->has_species()) {
          bcs().fillpatch(lev, time, BCFillVar::X, leveldata().X(), ng);
        }
        if (leveldata(lev)->has_tracer()) {
          bcs().fillpatch(lev, time, BCFillVar::tracer, leveldata().tracer(), ng);
        }
      }
      if (m_verbose > 0) { Print() << " Finished fillpatch before predictor\n"; }

      bool proj_2_corr = true;
      apply_corrector(stepData, m_timer.time(), m_timer.dt(), m_timer.prev_dt(), proj_2_corr);

      timings[t_fluid] += (ParallelDescriptor::second() - timer_corrector);
    }
  }


  ///////////////////////////////////////////////////////////////////////////////
  //                             ACTUAL PARTICLE STEP                          //
  ///////////////////////////////////////////////////////////////////////////////
  int nsubsteps(0);

  if (m_dem.solve() || m_pic.solve()) {

    EvolveSolids(m_timer.dt(), m_timer.time(), nsubsteps,
        timings[t_solids], timings[t_coupling], stepData.acceleration(),
        stepData.avg_particle_data_const(), stepData.eb_temperature());
  }




  ///////////////////////////////////////////////////////////////////////////////
  //                            REPORT STEP INFO                               //
  ///////////////////////////////////////////////////////////////////////////////
  ParallelDescriptor::ReduceRealMax(timings.dataPtr(), timings.size(),
    ParallelDescriptor::IOProcessorNumber());

  if (ParallelDescriptor::IOProcessor()) {
    if(fluid.solve())
      std::cout << "   Time per fluid step      " << timings[t_fluid] << '\n';

    if(m_dem.solve())
      std::cout << "   Time per " << nsubsteps
                << " particle steps " << timings[t_solids] << '\n';

    if(m_pic.solve())
      std::cout << "   Time per parcel step " << timings[t_solids] << '\n';

    if((m_dem.solve() || m_pic.solve()) && fluid.solve()) {
      std::cout << "   Coupling time per step   " << timings[t_coupling] << '\n';
    }
  }

  if (test_tracer_conservation) {

    // Compute the monitor on coarsest level
    const int lev = 0;
    EulerianMonitor::VolumeIntegral monitor(leveldata(), m_eb->factory(), geom[lev].ProbDomain());

    // we compute the monitor on the only component of the ep_g Multifab
    const int comp = 0;
    auto sum_trac_vol = monitor.volume_integral(lev, {leveldata().tracer(lev)}, {comp});
    auto sum_trac_weighted_vol = monitor.volume_weighted_integral(lev, {leveldata().tracer(lev)}, {comp});

    // sum_vol_orig is a vector with only one component
    const int var = 0;

    const Real* dx = geom[lev].CellSize();
    const Real cell_volume = dx[0] * dx[1] * dx[2];

    Print() << "Sum tracer volume wgt = "
            << sum_trac_vol[var] / cell_volume << " "
            << sum_trac_weighted_vol[var] / cell_volume << '\n';
  }
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  // TODO TODO TODO: check the following
//  std::swap(therm_p, therm_po);

  BL_PROFILE_REGION_STOP("mfix::Evolve");
}


// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix::
prepareToStep ( MFIXStepData& a_stepData )
{
  Real const time( m_timer.time() );
  int const avg_ncomp = avg_pc_parms().nComp();

  ///////////////////////////////////////////////////////////////////////////////
  //                          PREPARE FLUID FOR STEP                           //
  ///////////////////////////////////////////////////////////////////////////////
  if ( fluid.solve() ) {

    Print() << "\n ============   NEW TIME STEP   ============ \n";

    a_stepData.define( time, fluid.solve_density(), fluid.solve_enthalpy(),
        (fluid.solve_species() ? fluid.nspecies() : 0), fluid.ntracer(),
        m_embedded_boundaries.has_flow(),
        m_embedded_boundaries.has_temperature());

    if ( advection_type() == AdvectionType::MOL )
    { a_stepData.include_corrector(); }

    if (m_coupling.include_virtual_mass())
    { a_stepData.include_acceleration(); }

    a_stepData.include_avg_particle_data(avg_ncomp);
  }


  ///////////////////////////////////////////////////////////////////////////////
  //                          PREPARE SOLIDS FOR STEP                          //
  ///////////////////////////////////////////////////////////////////////////////
  // Sort particles by cell to improve in-memory locality
  if ( has_particles() ) {

    pc->sortNow( m_timer.nstep() );

    if ( avg_ncomp > 0 ) {
      average_pc_data_to_grid( a_stepData.avg_particle_data() );
    }
  }


  ///////////////////////////////////////////////////////////////////////////////
  //                          Impose boundary conditions                       //
  ///////////////////////////////////////////////////////////////////////////////
  if ( fluid.solve() ) {

    for (int lev(0); lev < nlev(); ++lev) {

      // Fill ghost nodes and reimpose boundary conditions
      //bcs().set_velocity_bcs(lev, time, leveldata().vel(lev));

      bcs().set_density_bcs(lev, time, leveldata().rho(lev));

      if (leveldata(lev)->has_species()) {
        bcs().set_species_bcs(lev, time, fluid.nspecies(), leveldata().X(lev));
      }

      // TODO: commenting the following makes BENCH03 GPU to pass
      //if (leveldata(lev)->has_tracer()) {
      //  bcs().set_tracer_bcs(lev, time, fluid, leveldata().tracer(lev));
      //}

      if (leveldata(lev)->has_temperature()) {
        bcs().set_temperature_bcs(lev, time, fluid, leveldata().T(lev));
      }

      if ( a_stepData.has_eb_temperature() ) {

        int const avgTp_comp = avg_pc_parms().comp(SoArealData::temperature);

        //const std::unique_ptr<MultiFab> avg_Tp = (comp < 0) ? nullptr :
        //    std::make_unique<MultiFab>(avg_data, amrex::make_alias, comp, 1);

        bcs().set_eb_temperature_bcs(lev, time, a_stepData.eb_temperature(),
            leveldata().epf_const(lev), leveldata().T_const(lev),
            avgTp_comp, a_stepData.avg_particle_data()[lev]);

      }
    }
  }

}
