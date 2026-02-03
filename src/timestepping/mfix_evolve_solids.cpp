#include <AMReX_BCUtil.H>

#include <mfix.H>
#include <mfix_fluid.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_reactions.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void mfix::
EvolveSolids ( Real const a_dt, Real const a_time, int& a_nsubsteps,
               Real& a_particle_timing, Real& a_coupling_timing,
               Vector< MultiFab* > const& a_DufDt,
               Vector< amrex::MultiFab const*> const& a_avg_particle_data,
               Vector< amrex::MultiFab      *> const& a_T_eb)
{
  BL_PROFILE_REGION_START("mfix::EvolveSolids");

  // Compute coupling forces with fluid for partlce updates
  if ( fluid.solve() ) {

    Real const timer_coupling( ParallelDescriptor::second() );

    Vector< MultiFab* > divtau(nlev(), nullptr);
    Vector< MultiFab* > eb_flow_vel(nlev(), nullptr);

    if (m_coupling.include_divtau()) {

      for (int lev(0); lev<nlev(); ++lev) {
        divtau[lev] = new MultiFab(grids[lev], dmap[lev],
          AMREX_SPACEDIM, 1, MFInfo(), *(m_eb->factory()[lev]));
        divtau[lev]->setVal(0.);

        if (m_embedded_boundaries.has_flow()) {

          eb_flow_vel[lev] = new MultiFab(grids[lev], dmap[lev],
            AMREX_SPACEDIM, nghost_state(), MFInfo(), *(m_eb->factory()[lev]));
          eb_flow_vel[lev]->setVal(0.0);

          bcs().set_eb_velocity_bcs(a_time, eb_flow_vel);
        }
      } // lev

      int const include_eddy_viscosity(1);
      diffOpVel()->computeDivTau(divtau, leveldata().vel(),
        leveldata().epf_const(), leveldata().rho_const(), leveldata().T_const(),
        a_avg_particle_data, leveldata().X_const(),
        GetVecOfConstPtrs(eb_flow_vel), include_eddy_viscosity);

      int const foextrap_bcs[AMREX_SPACEDIM] =
        {BCType::foextrap, BCType::foextrap, BCType::foextrap};

      Vector<BCRec> divtau_bcs(AMREX_SPACEDIM);
      divtau_bcs[0] = BCRec(foextrap_bcs, foextrap_bcs);
      divtau_bcs[1] = BCRec(foextrap_bcs, foextrap_bcs);
      divtau_bcs[2] = BCRec(foextrap_bcs, foextrap_bcs);

      for (int lev(0); lev<nlev(); ++lev) {

        divtau[lev]->FillBoundary(geom[lev].periodicity());

        FillDomainBoundary(*divtau[lev], geom[lev], divtau_bcs);

        EB_set_covered(*divtau[lev], 0, AMREX_SPACEDIM, 1, covered_val);
      }
    } // drag includes divtau

    // This returns the drag force on the particle
    Real new_time = a_time+a_dt;

    calc_txfr_particle(new_time, a_dt, leveldata().vel(), leveldata().T(),
      leveldata().grad_p_const(), GetVecOfConstPtrs(divtau),
      GetVecOfConstPtrs(a_DufDt));

    // Clean up local MultiFab pointers
    for (int lev(0); lev < nlev(); ++lev) {
      if (divtau[lev] != nullptr) { delete divtau[lev]; }
      if (eb_flow_vel[lev] != nullptr) { delete eb_flow_vel[lev]; }
    }

    a_coupling_timing += (ParallelDescriptor::second() - timer_coupling);

  }




  /****************************************************************************
   *                                                                          *
   * Evolve Particles (Using Particle MD)                                     *
   *                                                                          *
   * The following logic is a little confusing and I'm not sure why we go to  *
   * the trouble of having two nearly identical calls into EvolveParticles.   *
   *                                                                          *
   * 1. I think we should only generate one levelset at the finest level      *
   *    specified by the user, and this should be passed to on all levels.    *
   *                                                                          *
   * 2. The second loop--while it may be needed if we ever start pushing      *
   *    particles to higher levels--should be over the number of levels in    *
   *    the particle container and NOT the number of grid levels.             *
   *                                                                          *
   ***************************************************************************/

  Real const timer_solids( ParallelDescriptor::second() );

  if (m_dem.solve()) {

    BL_PROFILE_REGION("DEM PARTICLE SOLVE");

    if (nlev() == 1) {

      int ilev = 0;

      MultiFab* T_eb = a_T_eb[ilev];

      if ( eb_parms().has_temperature() ) {

        AMREX_ASSERT( a_T_eb[ilev] != nullptr );

        const DistributionMapping& Teb_dm = a_T_eb[ilev]->DistributionMap();
        const BoxArray&            Teb_ba = a_T_eb[ilev]->boxArray();

        const DistributionMapping& pc_dm = pc->ParticleDistributionMap(ilev);
        const BoxArray&            pc_ba = pc->ParticleBoxArray(ilev);

        int ncomp = a_T_eb[ilev]->nComp();
        int ngrow = 1;

        // Different grids, we need to copy from a_T_eb into T_eb
        if ((Teb_dm != pc_dm) || (Teb_ba != pc_ba)) {

          T_eb = new MultiFab(Teb_ba, pc_dm, ncomp, ngrow, MFInfo(),
              *(eb()->particle_factory()[ilev]) );

          T_eb->ParallelCopy(*a_T_eb[ilev], 0, 0, ncomp, 1, 1);

        }
      }

      //___________________________________________________________________
      // Single level case: the refined level-set is stored on level 1,
      // everything else lives on level 0

      const MultiFab* ls_data = m_eb->level_sets()[1].get();

      pc->EvolveParticles(ilev, a_dt, mfix::gravity,
            m_eb->factory()[ilev].get(), m_eb->particle_factory()[ilev].get(),
            ls_data, m_eb->levelset_refinement(), m_loadbalance.get(),
            T_eb, a_nsubsteps);

      // Clean up the local copy of T_eb
      if (T_eb != a_T_eb[ilev]) { delete T_eb; }

    } else {

      //___________________________________________________________________
      // Multi-level case: each level is treated separately

      // This loop should be over particle container levels
      for (int lev(0); lev<nlev(); ++lev) {

        MultiFab* T_eb = a_T_eb[lev];

        if ( eb_parms().has_temperature() ) {

          const DistributionMapping& Teb_dm = a_T_eb[lev]->DistributionMap();
          const BoxArray&            Teb_ba = a_T_eb[lev]->boxArray();

          const DistributionMapping& pc_dm = pc->ParticleDistributionMap(lev);
          const BoxArray&            pc_ba = pc->ParticleBoxArray(lev);

          int ncomp = a_T_eb[lev]->nComp();
          int ngrow = 1;

          // Different grids, we need to copy from a_T_eb into T_eb
          if ((Teb_dm != pc_dm) || (Teb_ba != pc_ba)) {

            T_eb = new MultiFab(Teb_ba, pc_dm, ncomp, ngrow, MFInfo(),
                *(eb()->particle_factory()[lev]) );

            T_eb->ParallelCopy(*a_T_eb[lev], 0, 0, ncomp, 1, 1);

          }
        }
        const MultiFab* ls_data = m_eb->level_sets()[lev].get();

        pc->EvolveParticles(lev, a_dt, mfix::gravity,
              m_eb->factory()[lev].get(), m_eb->particle_factory()[lev].get(),
              ls_data, 1, m_loadbalance.get(), a_T_eb[lev], a_nsubsteps);

        // Clean up the local copy of T_eb
        if (T_eb != a_T_eb[lev]) { delete T_eb; }

      }
    }

    // This loop should be over particle container levels
    if (m_rw->report_mass_balance) {
      for (int lev(0); lev < nlev(); lev ++ ) {
        m_mass_balance->ComputeMassOutflow(lev);
      }
    }

    // This loop should be over particle container levels
    for (int lev(0); lev < nlev(); lev ++ ) {
      pc->PostEvolveParticles(lev);
    }

  } // end DEM solve

  if (m_pic.solve()) {

      BL_PROFILE_REGION("PIC PARTICLE SOLVE");
      //const IntVect min_epg_cell = m_rw->mfix_print_min_epg();
      EvolveParcels(a_dt, a_time, mfix::gravity, m_eb->levelset_refinement(),
          m_loadbalance.get());

      if (m_rw->report_mass_balance) {
        // This loop should be over particle container levels
        for (int lev(0); lev < nlev(); lev ++ ) {
          m_mass_balance->ComputeMassOutflow(lev);
        }
      }

      PostEvolveParcels();
  }


  if (m_dem.solve() || m_pic.solve()) {
    if (pc->UseConstraint()) {
      // This loop should be over particle container levels
      for (int lev = 0; lev < nlev(); lev ++ ) {
        pc->ImposeMean(lev);
      }
    }
  }

  a_particle_timing = (ParallelDescriptor::second() - timer_solids);
}
