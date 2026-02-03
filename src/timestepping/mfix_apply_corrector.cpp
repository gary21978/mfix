#include <AMReX_VisMF.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <mfix.H>
#include <mfix_run_on.H>
#include <mfix_mf_helpers.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_solvers.H>

#include <mfix_fluid_update.H>

using namespace Solvers;

/***************************************************************************
 * MOL corrector step                                                      *
 *                                                                         *
 *  1. Compute the MAC projected velocities, (εu)_mac                      *
 *                                                                         *
 *  2. Compute the the convective terms:                                   *
 *                                                                         *
 *     density:   A_ρ  = ∇·(εu)_mac ρ                                      *
 *     species:   A_Xk = ∇·(εu)_mac (ρ Xk)                                 *
 *     enthalpy:  A_h  = ∇·(εu)_mac (ρ h)                                  *
 *     tracer:    A_s  = ∇·(εu)_mac (ρ s)                                  *
 *                                                                         *
 *     velocity:  A_u  = (εu)_mac · ∇ u    -- convective differencing      *
 *                A_ρu = ∇·(εu)_mac (ρ u)  -- advecting momentum           *
 *                                                                         *
 *  3. Compute updates: density, species, enthalpy, temperature, tracers   *
 *     and velocity.                                                       *
 *                                                                         *
 *     ρ^(n+1), Xk^(n+1), h^(n+1), T^(n+1), s^(n+1), u^**                  *
 *                                                                         *
 *  4. Apply nodal projection:  ∇·(ε/ρ)∇φ =  ∇·(εu^**) - ∇·(εu)            *
 *                                                                         *
 *     u^(n+1) = u^** - (1/ρ)∇φ                                            *
 *     pert_p = (1/Δt) ∇φ                                                  *
 *                                                                         *
 ***************************************************************************/
void mfix::
apply_corrector ( MFIXStepData& a_stepData, Real time,
                  Real l_dt, Real l_prev_dt, bool proj_2)
{
    BL_PROFILE("mfix::apply_corrector");

    StepData& predictor = a_stepData.predictor();
    StepData& corrector = a_stepData.corrector();

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + l_dt;

    InterphaseChemTxfrIndexes chem_txfr_idxs(fluid.nspecies(), reactions.solve());

    const bool explicit_diffusive_enthalpy = false;
    const bool explicit_diffusive_trac     = false;
    const bool explicit_diffusive_species  = false;


    // *************************************************************************************
    // Compute explicit diffusive terms
    // *************************************************************************************
    {
      const bool constraint = !(fluid.constraint.isIncompressibleFluid());

      for ( int lev(0); lev < nlev(); ++lev) {

        bcs().set_density_bcs(lev, time, leveldata().rho(lev));

        if (leveldata(lev)->has_enthalpy()) {
          bcs().set_enthalpy_bcs(lev, time, fluid, leveldata().h(lev));
        }
        if (leveldata(lev)->has_temperature()) {
          bcs().set_temperature_bcs(lev, time, fluid, leveldata().T(lev));
        }
        if (leveldata(lev)->has_species()) {
          bcs().set_species_bcs(lev, time, fluid.nspecies(), leveldata().X(lev));
        }
      }

      if (fluid.solve_species() && (explicit_diffusive_species || constraint)) {
        diffOpSpecies()->computeFlux(a_stepData.J_k(), leveldata().X(),
            leveldata().epf_const(), leveldata().rho_const());

        diffOpSpecies()->computeDivJ(corrector.div_J(), a_stepData.J_k());
      }

      if (fluid.solve_enthalpy() && (explicit_diffusive_enthalpy || constraint)) {
        diffOpEnergy()->computeLapT(corrector.lap_T(), leveldata().T(),
            leveldata().epf_const(), a_stepData.eb_temperature_const());
      }

      if (fluid.solve_enthalpy() && fluid.solve_species()) {

        diffOpEnergy()->computeDivhJ(corrector.div_hJ(),
            a_stepData.h_k(), a_stepData.J_k(), leveldata().T_const(),
            /*update_enthalpies=*/ 1);
      }

      if (fluid.solve_tracer() && explicit_diffusive_trac) {
        diffOpTracer()->computeLap(corrector.lap_tracer(), leveldata().tracer(),
            leveldata().epf_const(), leveldata().rho_const());
      }

      // Call the bc routines again to enforce the ext_dir condition on the faces
      // (the diffusion operator can move those to ghost cell centers)
      for ( int lev(0); lev < nlev(); ++lev) {

        if (leveldata(lev)->has_enthalpy()) {
          bcs().set_enthalpy_bcs(lev, time, fluid, leveldata().h(lev));
        }
        if (leveldata(lev)->has_temperature()) {
          bcs().set_temperature_bcs(lev, time, fluid, leveldata().T(lev));
        }
        if (leveldata(lev)->has_species()) {
          bcs().set_species_bcs(lev, time, fluid.nspecies(), leveldata().X(lev));
        }
      }
    }

    // *************************************************************************************
    // Compute RHS for the MAC projection
    // *************************************************************************************
    compute_MAC_proj_RHS( a_stepData.RHS_proj(), a_stepData.depdt_const(),
        leveldata().epf_const(), leveldata().rho_const(),
        leveldata().X_const(), leveldata().T_const(),
        corrector.lap_T_const(), corrector.div_hJ_const(), corrector.div_J_const(),
        leveldata().txfr(), leveldata().chem_txfr_const(),
        therm_p, corrector.RHS_therm_p);

    // *************************************************************************************
    // Compute the explicit advective terms
    // *************************************************************************************
    compute_MAC_projected_velocities(time, l_dt, leveldata().vel_const(),
        a_stepData.RHS_proj_const(), a_stepData.mac_vel(), leveldata().mac_phi(),
        leveldata().epf_const(), leveldata().rho_const(),
        leveldata().txfr_const(), a_stepData.eb_flow_velocity_const(),
        a_stepData.vel_forces(), leveldata().div_tau_const());

    compute_convective_term(corrector.conv_velocity(), corrector.conv_scalars(),
        corrector.conv_species(), a_stepData.vel_forces(), a_stepData.tracer_forces(),
        leveldata().vel_const(), leveldata().epf_const(), leveldata().rho_const(),
        leveldata().h_const(), leveldata().tracer_const(), leveldata().X_const(),
        leveldata().txfr_const(), a_stepData.eb_flow_velocity_const(),
        a_stepData.eb_flow_scalars_const(), a_stepData.eb_flow_species_const(),
        a_stepData.ep_u_mac(), a_stepData.ep_v_mac(), a_stepData.ep_w_mac(), l_dt, new_time);

    FluidUpdate updateFluid( nlev(), StepType::Corrector, geom, grids, fluid,
        leveldata(), a_stepData, bcs(), reactions.solve(), time, l_dt);

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    updateFluid.Density( );

    // *************************************************************************
    // Update species mass fraction
    // *************************************************************************
    if (fluid.solve_species()) {

      updateFluid.Species( diffOpSpecies() );

    } // fluid.solve_species


    if (fluid.solve_enthalpy()) {

      updateFluid.Energy( therm_po, therm_p, diffOpEnergy(),
          newton_abstol, newton_reltol, newton_maxiter);

      // ***********************************************************************
      // Add the drag and enthalpy terms implicitly
      // ***********************************************************************
      if (m_dem.solve() || m_pic.solve())
        add_enthalpy_txfr_implicit(l_dt, leveldata().h(), leveldata().T(), leveldata().X_const(),
            leveldata().txfr_const(), leveldata().rho_const(), leveldata().epf_const());

    } // fluid.solve_enthalpy


    // *************************************************************************************
    // Update tracer(s)
    // *************************************************************************************
    if (fluid.solve_tracer()) { updateFluid.Tracers( fluid.ntracer(), diffOpTracer() ); }


    // *************************************************************************************
    // Define (or if advection_type != "MOL", re-define) the forcing terms, without the
    //    viscous terms and using the half-time density
    // *************************************************************************************
    compute_vel_forces(a_stepData.vel_forces(), leveldata().vel_const(),
        a_stepData.rho_nph_const(), leveldata().txfr_const());

    updateFluid.ComputeVelTxfrSrc(advection_type(), m_advect_momentum,
        m_coupling.include_divtau());

    m_porous_media.ComputeVelSource( a_stepData.vel_S_p(), a_stepData.vel_S_c(),
        leveldata().rho_const(), leveldata().vel_const(),
        leveldata().T_const(), leveldata().X_const(),
        fluid.props.data<run_on>(), /*S_fac = */ 0.5);

    // *************************************************************************************
    // Update velocity with convective update, diffusive update, gp and gravity source terms
    // *************************************************************************************
    int const include_src( m_dem.solve() || m_pic.solve() || (m_porous_media.nregions() > 0));

    updateFluid.Velocity( m_advect_momentum, include_src,
        m_coupling.include_divtau(), diffOpVel());

    // *************************************************************************************
    // Apply projection
    // *************************************************************************************
    Real RHS_therm_p = 0.5*(predictor.RHS_therm_p +
                            corrector.RHS_therm_p);

    compute_nodal_proj_RHS(a_stepData.RHS_proj(), a_stepData.depdt_const(),
        leveldata().epf_old_const(), leveldata().epf_const(),
        leveldata().rho_old_const(), leveldata().rho_const(),
        leveldata().h_old_const(), leveldata().h_const(),
        leveldata().T_const(), corrector.conv_scalars_const(),
        leveldata().X_old_const(), leveldata().X_const(), corrector.conv_species_const(),
        therm_p, RHS_therm_p, m_timer.dt());

    apply_nodal_projection(new_time, l_dt, l_prev_dt, proj_2,
        a_stepData.RHS_proj(), leveldata().vel_old(), leveldata().vel(),
        leveldata().pert_p(), leveldata().grad_p(), leveldata().epf_const(),
        a_stepData.rho_nph_const(), a_stepData.eb_flow_velocity_const());


    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells( leveldata().vel(), a_stepData.ep_u_mac_const(),
         a_stepData.ep_v_mac_const(), a_stepData.ep_w_mac_const(),
         a_stepData.eb_flow_velocity_const());


    // *************************************************************************************
    // Compute fluid acceleration
    // *************************************************************************************
    if (m_coupling.include_virtual_mass()) {

      updateFluid.ComputeAcceleration( advection_type(), m_advect_momentum);
    }
}
