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
 * MOL predictor step // Godunov step:                                     *
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
 *     MOL:      ρ^*, Xk^*, h^*, T^*, s^*, u^**                            *
 *                                                                         *
 *     Godunov:  ρ^(n+1), Xk^(n+1), h^(n+1), T^(n+1), s^(n+1), u^**        *
 *                                                                         *
 *  4. Apply nodal projection:  ∇·(ε/ρ)∇φ =  ∇·(εu^*) - ∇·(εu)             *
 *                                                                         *
 *     u^(n+1) = u^** - (1/ρ)∇φ                                            *
 *     pert_p = (1/Δt) ∇φ                                                  *
 *                                                                         *
 ***************************************************************************/
void mfix::
apply_predictor ( MFIXStepData& a_stepData,
                  Real time, Real l_dt, Real l_prev_dt, bool proj_2)
{
    BL_PROFILE("mfix::apply_predictor");

    StepData& predictor = a_stepData.predictor();

    // Local flag for explicit diffusion
    bool const l_explicit_diff(predictor_diff_type() == DiffusionType::Explicit);

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + l_dt;

    InterphaseChemTxfrIndexes chem_txfr_idxs(fluid.nspecies(), reactions.solve());

    // *************************************************************************************
    // Compute explicit diffusive terms
    // *************************************************************************************
    {
      const bool constraint = !(fluid.constraint.isIncompressibleFluid());

      for ( int lev(0); lev < nlev(); ++lev) {

        bcs().set_density_bcs(lev, time, leveldata().rho_old(lev));

        if (leveldata(lev)->has_enthalpy()) {
          bcs().set_enthalpy_bcs(lev, time, fluid, leveldata().h_old(lev));
        }
        if (leveldata(lev)->has_temperature()) {
          bcs().set_temperature_bcs(lev, time, fluid, leveldata().T_old(lev));
        }
        if (leveldata(lev)->has_species()) {
          bcs().set_species_bcs(lev, time, fluid.nspecies(), leveldata().X_old(lev));
        }
      }

      if (need_divtau()) {
        diffOpVel()->computeDivTau(leveldata().div_tau(), leveldata().vel_old(), leveldata().epf_const(),
          leveldata().rho_old_const(), leveldata().T_old_const(), a_stepData.avg_particle_data_const(),
          leveldata().X_old_const(), a_stepData.eb_flow_velocity_const(), /*include_eddy_viscosity = */ 1);
      }


      if (fluid.solve_species()) {
        diffOpSpecies()->computeFlux(a_stepData.J_k(), leveldata().X_old(),
            leveldata().epf_const(), leveldata().rho_old_const());

        diffOpSpecies()->computeDivJ(predictor.div_J(),
            a_stepData.J_k());
      }


      if (fluid.solve_enthalpy() && (l_explicit_diff || constraint)) {
        diffOpEnergy()->computeLapT(predictor.lap_T(), leveldata().T_old(),
            leveldata().epf_const(), a_stepData.eb_temperature_const());
      }

      if (fluid.solve_enthalpy() && fluid.solve_species()) {

        diffOpEnergy()->computeDivhJ(predictor.div_hJ(),
            a_stepData.h_k(), a_stepData.J_k(), leveldata().T_old_const(),
            /*update_enthalpies=*/ 1);
      }

      if (fluid.solve_tracer() &&  l_explicit_diff) {
        diffOpTracer()->computeLap(predictor.lap_tracer(),
            leveldata().tracer_old(), leveldata().epf_const(), leveldata().rho_old_const());
      }

      // Call the bc routines again to enforce the ext_dir condition on the faces
      // (the diffusion operator may move those to ghost cell centers)
      for ( int lev(0); lev < nlev(); ++lev) {

        if (leveldata(lev)->has_enthalpy()) {
          bcs().set_enthalpy_bcs(lev, time, fluid, leveldata().h_old(lev));
        }
        if (leveldata(lev)->has_temperature()) {
          bcs().set_temperature_bcs(lev, time, fluid, leveldata().T_old(lev));
        }
        if (leveldata(lev)->has_species()) {
          bcs().set_species_bcs(lev, time, fluid.nspecies(), leveldata().X_old(lev));
        }
      }
    }

    // *************************************************************************************
    // Compute RHS for the MAC projection
    // *************************************************************************************
    compute_MAC_proj_RHS( a_stepData.RHS_proj(), a_stepData.depdt_const(),
        leveldata().epf_const(), leveldata().rho_old_const(),
        leveldata().X_old_const(), leveldata().T_old_const(),
        predictor.lap_T_const(), predictor.div_hJ_const(), predictor.div_J_const(),
        leveldata().txfr(), leveldata().chem_txfr_const(),
        therm_po, predictor.RHS_therm_p);


    // *************************************************************************************
    // Compute the explicit advective terms
    // Note that "conv_velocity" returns update to (epf u)
    // Note that "conv_scalar" returns update to (epf rho), (epf rho h) and (epf rho tracer)
    // *************************************************************************************
    compute_MAC_projected_velocities( time, l_dt, leveldata().vel_old_const(),
        a_stepData.RHS_proj_const(), a_stepData.mac_vel(), leveldata().mac_phi(),
        a_stepData.epf_nph_const(), leveldata().rho_old_const(),
        leveldata().txfr_const(), a_stepData.eb_flow_velocity_const(),
        a_stepData.vel_forces(), leveldata().div_tau_const());

    compute_convective_term(predictor.conv_velocity(), predictor.conv_scalars(),
        predictor.conv_species(), a_stepData.vel_forces(), a_stepData.tracer_forces(),
        leveldata().vel_old_const(), leveldata().epf_const(), leveldata().rho_old_const(),
        leveldata().h_old_const(), leveldata().tracer_old_const(), leveldata().X_old_const(),
        leveldata().txfr_const(), a_stepData.eb_flow_velocity_const(),
        a_stepData.eb_flow_scalars_const(), a_stepData.eb_flow_species_const(),
        a_stepData.ep_u_mac(), a_stepData.ep_v_mac(), a_stepData.ep_w_mac(), l_dt, time);

    FluidUpdate updateFluid( nlev(), StepType::Predictor, geom, grids, fluid,
        leveldata(), a_stepData, bcs(), reactions.solve(), time, l_dt);

    if (predictor_diff_type() == DiffusionType::Explicit)
    { updateFluid.setExplicitDiffusion(); }


    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    updateFluid.Density();


    // *************************************************************************
    // Update species mass fractions
    // *************************************************************************
    if (fluid.solve_species()) {

      updateFluid.Species( diffOpSpecies() );

      // Compute divJ for Crank-Nicolson update for species diffusion in corrector
      if (!l_explicit_diff) {
        diffOpSpecies()->computeDivJ(predictor.div_J(), a_stepData.J_k());
      }

    } // solve_species


    if (fluid.solve_enthalpy()) {

      updateFluid.Energy( therm_po, therm_p, diffOpEnergy(),
          newton_abstol, newton_reltol, newton_maxiter);

      // ***********************************************************************
      // Add the dconvective heat transfer terms implicitly to h
      // ***********************************************************************
      if (m_dem.solve() || m_pic.solve()) {
        add_enthalpy_txfr_implicit(l_dt, leveldata().h(), leveldata().T(), leveldata().X_const(),
            leveldata().txfr_const(), leveldata().rho_const(), leveldata().epf_const());
      }
    } // fluid.solve_enthalpy


    // *************************************************************************************
    // Update tracer(s)
    // *************************************************************************************
    if (fluid.solve_tracer()) {
      updateFluid.Tracers( fluid.ntracer(), diffOpTracer() );
    }


    // *************************************************************************************
    // Define (or if advection_type != "MOL", re-define) the forcing terms, without the
    //    viscous terms and using the half-time density
    // *************************************************************************************
    compute_vel_forces(a_stepData.vel_forces(), leveldata().vel_old_const(),
        a_stepData.rho_nph_const(), leveldata().txfr_const());

    updateFluid.ComputeVelTxfrSrc( advection_type(), m_advect_momentum,
        m_coupling.include_divtau());

    m_porous_media.ComputeVelSource( a_stepData.vel_S_p(), a_stepData.vel_S_c(),
        leveldata().rho_const(), leveldata().vel_const(),
        leveldata().T_const(), leveldata().X_const(),
        fluid.props.data<run_on>(), /*S_fac = */ 1.0);


    // *************************************************************************************
    // Update velocity with convective update, diffusive update, gp and gravity source terms
    // *************************************************************************************
    int const include_src( m_dem.solve() || m_pic.solve() || (m_porous_media.nregions() > 0));

    updateFluid.Velocity( m_advect_momentum, include_src,
        m_coupling.include_divtau(), diffOpVel());


    // *************************************************************************************
    // Project velocity field
    //
    //  1. Compute nodal projection RHS
    //  2. Apply nodal projectoin
    //
    // *************************************************************************************
    Real RHS_therm_p = predictor.RHS_therm_p;

    compute_nodal_proj_RHS(a_stepData.RHS_proj(), a_stepData.depdt_const(),
        leveldata().epf_old_const(), leveldata().epf_const(),
        leveldata().rho_old_const(), leveldata().rho_const(),
        leveldata().h_old_const(), leveldata().h_const(),
        leveldata().T_const(), predictor.conv_scalars_const(),
        leveldata().X_old_const(), leveldata().X_const(), predictor.conv_species_const(),
        therm_p, RHS_therm_p, m_timer.dt());

    apply_nodal_projection(new_time, l_dt, l_prev_dt, proj_2,
        a_stepData.RHS_proj(), leveldata().vel_old(), leveldata().vel(),
        leveldata().pert_p(), leveldata().grad_p(), leveldata().epf_const(),
        a_stepData.rho_nph_const(), a_stepData.eb_flow_velocity_const());


    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells (leveldata().vel(), a_stepData.ep_u_mac_const(),
         a_stepData.ep_v_mac_const(), a_stepData.ep_w_mac_const(),
         a_stepData.eb_flow_velocity_const());

    // *************************************************************************************
    // Compute fluid acceleration
    // *************************************************************************************
    if (m_coupling.include_virtual_mass()) {
      updateFluid.ComputeAcceleration(advection_type(), m_advect_momentum);
    }
}
