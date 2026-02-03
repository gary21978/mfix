#include <AMReX_BCUtil.H>

#include <mfix.H>
#include <mfix_utils.H>

void mfix::
Evolve_MOL ()
{
  BL_PROFILE_REGION_START("mfix::Evolve_MOL");

  int const t_fluid(0), t_solids(1), t_coupling(2);
  Vector<Real> timings = {0., 0., 0.};

  ///////////////////////////////////////////////////////////////////////////////
  //                           PREPARE FLUID FOR STEP                          //
  ///////////////////////////////////////////////////////////////////////////////
  MFIXStepData stepData(nlev(), grids, dmap, eb()->factory_const(),
      nghost_state(), nghost_force(), m_boundary_conditions );

  prepareToStep( stepData );

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

    if (m_dem.solve() || m_pic.solve())
    { m_coupling.calc_transfer_coeffs(fluid, pc, geom[0].data()); }

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
  //                    DEPOSIT FLUID-SOLIDS COUPLING TO GRID                  //
  ///////////////////////////////////////////////////////////////////////////////
  if (fluid.solve() && (m_dem.solve() || m_pic.solve()) )
  { deposit_forces_to_grid(m_timer.dt()); }


  ///////////////////////////////////////////////////////////////////////////////
  //                        EXTRAPOLATE VOLUME FRACTION                        //
  ///////////////////////////////////////////////////////////////////////////////
  if ( fluid.solve() ) {

    Real const timer_volfrac( ParallelDescriptor::second() );

    if (m_dem.solve() || m_pic.solve()) {

      for (int lev(0); lev < nlev(); ++lev) {

        // We haven't computed the new dt yet, so timer.dt() contains the
        // previous time step value: Δt^(n-1)
        if ( !almostEqual(m_timer.dt(), Real(0.)) ) {

          MultiFab::LinComb(*stepData.depdt(lev), 1.0, *leveldata().epf(lev), 0,
              -1.0, *leveldata().epf_old(lev), 0, 0, 1, stepData.depdt(lev)->nGrow());

          Real const inv_dt( 1_rt/m_timer.dt() );
          stepData.depdt(lev)->mult(inv_dt, 0, 1, stepData.depdt(lev)->nGrow());
        }

        MultiFab::Copy(*leveldata().epf_old(lev), *leveldata().epf(lev),
            0, 0, 1, leveldata().epf(lev)->nGrow());

        MultiFab::Copy(*stepData.epf_nph(lev), *leveldata().epf(lev),
            0, 0, 1, stepData.epf_nph(lev)->nGrow());
      }

      //const IntVect min_epg_cell = m_rw->mfix_print_min_epg();
      timings[t_coupling] = (ParallelDescriptor::second() - timer_volfrac);

    } else { // no dem or pic

      if (m_porous_media.nregions() > 0) {
        m_porous_media.impose_vfrac(leveldata().epf());
      }
    }
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

    std::swap(therm_p, therm_po);

    apply_predictor(stepData, m_timer.time(), m_timer.dt(),
        m_timer.prev_dt(), /*proj_2_pred=*/true);

    timings[t_fluid] += (ParallelDescriptor::second() - timer_predictor);
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
  //                           COMPUTE VOLUME FRACTION                         //
  ///////////////////////////////////////////////////////////////////////////////

  if ( fluid.solve() ) {

    Real const timer_volfrac( ParallelDescriptor::second() );

    if (m_dem.solve() || m_pic.solve()) {

      deposit_volume_to_grid( m_timer.time() );

      if (m_porous_media.nregions() > 0) {
        m_porous_media.impose_vfrac(leveldata().epf());
      }

      for (int lev(0); lev < nlev(); ++lev) {

        // Set the volume fraction to the correct values in inflow faces.
        bcs().set_epf_bcs(lev, m_timer.time(), leveldata().epf(lev), BCType::ext_dir);

        // Compute the 'half-time' volume volume fraction as an average of old
        // and new fields. The old value was used for the predictor so here we
        // add the new field and divide by two.

        MultiFab::Add(*stepData.epf_nph(lev), *leveldata().epf(lev),
            0, 0, 1, leveldata().epf(lev)->nGrow());

        stepData.epf_nph(lev)->mult(0.5, 0, 1, leveldata().epf(lev)->nGrow());

        if ( !almostEqual(m_timer.dt(), Real(0.)) ) {

          int const ngrow( stepData.depdt(lev)->nGrow() );

          // depdt contains (ε^n - ε^(n-1)) / dt^(n-1)

          MultiFab depf_new( grids[lev], dmap[lev], /*ncomp=*/1, ngrow,
              MFInfo(), m_eb->Factory(lev) );

          // depf_new = (ε^(n+1) - ε^n)
          MultiFab::LinComb( depf_new, 1.0, *leveldata().epf(lev), 0,
              -1.0, *leveldata().epf_old(lev), 0, 0, 1, ngrow);

          Real const coeff( -Real(3.0)/m_timer.dt() );

          // depdt = depdt - 3(depf)/dt^n
          // = (ε^n - ε^(n-1)) / dt^(n-1) - 3(ε^(n+1) - ε^n) / dt^n
          MultiFab::Saxpy( *stepData.depdt(lev), coeff, depf_new,
            /*scomp*/0, /*dcomp*/0, /*ncomp*/1, ngrow);

          // depdt = (1/2)( 3(ε^(n+1) - ε^n) / dt^n - (ε^n - ε^(n-1)) / dt^(n-1))
          stepData.depdt(lev)->mult(Real(-0.5), 0, 1, ngrow);
        }

      }
      //const IntVect min_epg_cell = m_rw->mfix_print_min_epg();
      timings[t_coupling] = (ParallelDescriptor::second() - timer_volfrac);

    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  //                    DEPOSIT FLUID-SOLIDS COUPLING TO GRID                  //
  ///////////////////////////////////////////////////////////////////////////////
  if (fluid.solve() && (m_dem.solve() || m_pic.solve()) )
  { deposit_forces_to_grid(m_timer.dt()); }

  ///////////////////////////////////////////////////////////////////////////
  //                           CORRECTOR STEP                              //
  ///////////////////////////////////////////////////////////////////////////
  if (fluid.solve()) {

    Real const timer_corrector(ParallelDescriptor::second());

    //  Do fillpatch in here
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
    }

    bool proj_2_corr = true;
    apply_corrector(stepData, m_timer.time(), m_timer.dt(), m_timer.prev_dt(), proj_2_corr);

    timings[t_fluid] += (ParallelDescriptor::second() - timer_corrector);
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

    // we compute the monitor on the only component of the epf Multifab
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
