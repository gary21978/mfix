#include <AMReX_VisMF.H>

#include <mfix.H>

using namespace amrex;

void
mfix::
project_velocity ()
{
  // Project velocity field to make sure initial velocity is divergence-free
  Real dummy_dt = 1.0;

  amrex::Print() << "Initial projection:\n";

  bool proj_2 = (!run_type(RunType::PIC2DEM)) ? true : false;
  Real time = 0.0;

  // Apply projection -- RHS=0 for now
  Vector< std::unique_ptr<MultiFab>> eb_flow_vel(nlev());
  Vector< std::unique_ptr<MultiFab>> S_cc(nlev());

  for (int lev(0); lev < nlev(); ++lev) {

    S_cc[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev],
        /*ncomp=*/1, /*ngrow=*/1, MFInfo(), m_eb->Factory(lev));
    S_cc[lev]->setVal(0.);

    if (m_embedded_boundaries.has_flow()) {
      eb_flow_vel[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev],
          /*ncomp=*/3, nghost_state(), MFInfo(), m_eb->Factory(lev));
      eb_flow_vel[lev]->setVal(0.0);
    }

    leveldata().grad_p(lev)->setVal(0);
  }

  apply_nodal_projection(time, dummy_dt, dummy_dt, proj_2,
      GetVecOfPtrs(S_cc), leveldata().vel_old(), leveldata().vel(),
      leveldata().pert_p(), leveldata().grad_p(), leveldata().epf_const(),
      leveldata().rho_const(), GetVecOfConstPtrs(eb_flow_vel));

  // We initialize p_g and gp back to zero
  for (int lev = 0; lev < nlev(); lev++) {
    leveldata().pert_p(lev)->setVal(0);
    leveldata().grad_p(lev)->setVal(0);
  }

  // Call the initial redistribution after the initial projection.
  if (m_redistribute_after_initial_nodal_proj) {
    InitialRedistribution(time);
  }
}



void mfix::initial_iterations ()
{
  MFIXStepData stepData(nlev(), grids, dmap, eb()->factory_const(),
      nghost_state(), nghost_force(), m_boundary_conditions);

  stepData.define( 0.0, fluid.solve_density(), fluid.solve_enthalpy(),
      (fluid.solve_species() ? fluid.nspecies() : 0), fluid.ntracer(),
      m_embedded_boundaries.has_flow(),
      m_embedded_boundaries.has_temperature());

  if (m_coupling.include_virtual_mass())
  { stepData.include_acceleration(); }

  int const avg_ncomp = avg_pc_parms().nComp();
  stepData.include_avg_particle_data(avg_ncomp);

  if ( has_particles() && (avg_ncomp > 0) ) {
    average_pc_data_to_grid( stepData.avg_particle_data() );
  }

  compute_dt(stepData.avg_particle_data_const());

  Print() << "Doing initial pressure iterations with dt = " << m_timer.dt() << "\n";

  // Fill ghost nodes and reimpose boundary conditions
  for (int lev(0); lev < nlev(); lev++) {

    bcs().set_velocity_bcs(lev, m_timer.time(), leveldata().vel(lev));
    bcs().set_density_bcs(lev, m_timer.time(), leveldata().rho(lev));

    if (leveldata(lev)->has_tracer()) {
      bcs().set_tracer_bcs(lev, m_timer.time(), fluid, leveldata().tracer(lev));
    }

    if (leveldata(lev)->has_temperature()) {
      bcs().set_temperature_bcs(lev, m_timer.time(), fluid, leveldata().T(lev));
    }

    if (leveldata(lev)->has_enthalpy()) {
      bcs().set_enthalpy_bcs(lev, m_timer.time(), fluid, leveldata().h(lev));
    }

    if (leveldata(lev)->has_species()) {
      bcs().set_species_bcs(lev, m_timer.time(), fluid.nspecies(), leveldata().X(lev));
    }

    if ( stepData.has_eb_temperature() ) {
      int const avgTp_comp = avg_pc_parms().comp(SoArealData::temperature);

      bcs().set_eb_temperature_bcs( lev, m_timer.time(), stepData.eb_temperature(),
          leveldata().epf_const(lev), leveldata().T_const(lev),
          avgTp_comp, stepData.avg_particle_data()[lev]);
    }
  }

  if (max_level > 0) {

    for (int lev(0); lev < nlev(); lev++) {

      Real const time = timer().time();

      bcs().fillpatch(lev, time, BCFillVar::vel, leveldata().vel(), 3);
      bcs().fillpatch(lev, time, BCFillVar::rho, leveldata().rho(), 3);

      if (leveldata(lev)->has_enthalpy()) {
        bcs().fillpatch(lev, time, BCFillVar::h, leveldata().h(), 3);
      }
      if (leveldata(lev)->has_temperature()) {
        bcs().fillpatch(lev, time, BCFillVar::T, leveldata().T(), 3);
      }
      if (leveldata(lev)->has_species()) {
        bcs().fillpatch(lev, time, BCFillVar::X, leveldata().X(), 3);
      }
      if (leveldata(lev)->has_tracer()) {
        bcs().fillpatch(lev, time, BCFillVar::tracer, leveldata().tracer(), 3);
      }
    }
  }

  // Copy vel into vel_old
  for (int lev = 0; lev < nlev(); lev++) {
    MultiFab::Copy(*(leveldata().vel_old(lev)), *(leveldata().vel(lev)), 0, 0,
        leveldata().vel(lev)->nComp(), leveldata().vel(lev)->nGrow());
  }

  ///////////////////////////////////////////////////////////////////////////////
  //                 COMPUTE FLUID-SOLIDS COUPLING AND REACTIONS               //
  ///////////////////////////////////////////////////////////////////////////////
  if (m_dem.solve() || m_pic.solve() ) {

    // Copy epf new into volfrac_nph
    for (int lev = 0; lev < nlev(); lev++) {

      MultiFab::Copy(*(leveldata().epf_old(lev)), *(leveldata().epf(lev)),
          0, 0, 1, leveldata().epf(lev)->nGrow());

      MultiFab::Copy(*stepData.epf_nph(lev), *(leveldata().epf(lev)),
            0, 0, 1, leveldata().epf(lev)->nGrow());

      // We copy the value inside the domain to the outside to avoid unphysical
      // volume fractions.
      bcs().set_epf_bcs(lev, m_timer.time(), leveldata().epf(lev), BCType::foextrap);
    }

    // Setup the interpolation MultiFab
    m_coupling.setup( Geom(), leveldata().epf(), leveldata().rho(), leveldata().vel(),
        leveldata().T(), leveldata().h(), leveldata().X(), eb()->pc_factory_const());

    // Reset the volume fractions back to the correct values at inflow faces.
    for (int lev = 0; lev < nlev(); lev++) {
      bcs().set_epf_bcs(lev, m_timer.time(), leveldata().epf(lev), BCType::ext_dir);
    }

    if (m_dem.solve() || m_pic.solve()) {
      m_coupling.calc_transfer_coeffs(fluid, pc, geom[0].data());
      deposit_forces_to_grid( m_timer.dt() );
    }

    // Not including reactions in initial iterations

    m_coupling.reset();
  }


  for (int iter(0); iter < m_initial_iterations; ++iter) {

    Print() << '\n' << "In initial_iterations: iter = " << iter << '\n';

    bool proj_2 = false;

    auto dt_copy = m_timer.dt();

    apply_predictor(stepData, m_timer.time(), m_timer.dt(), dt_copy, proj_2);

    therm_p = therm_po;
    for (int lev(0); lev<nlev(); lev++) {

      leveldata(lev)->resetNewWithOld();

      // Reset the boundary values (necessary if they are time-dependent)
      bcs().set_velocity_bcs(lev, m_timer.time(), leveldata().vel(lev));
      bcs().set_density_bcs(lev, m_timer.time(), leveldata().rho(lev));

      if (leveldata(lev)->has_tracer()) {
        bcs().set_tracer_bcs(lev, m_timer.time(), fluid, leveldata().tracer(lev));
      }
      if (leveldata(lev)->has_temperature()) {
        bcs().set_temperature_bcs(lev, m_timer.time(), fluid, leveldata().T(lev));
      }
      if (leveldata(lev)->has_enthalpy()) {
        bcs().set_enthalpy_bcs(lev, m_timer.time(), fluid, leveldata().h(lev));
      }
      if (fluid.solve_species()) {
        bcs().set_species_bcs(lev, m_timer.time(), fluid.nspecies(), leveldata().X(lev));
      }
    }

  } // initial iterations

}
