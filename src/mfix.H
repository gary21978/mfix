#ifndef MFIX_LEVEL_H_
#define MFIX_LEVEL_H_

#include <iostream>
#include <memory>
#include <string>


#include <AMReX_AmrCore.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2_IndexSpace_STL.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <hydro_NodalProjector.H>
#include <AMReX_PhysBCFunct.H>
#include <hydro_MacProjector.H>

#include <mfix_diffusion_op.H>
#include <mfix_pc.H>
#include <mfix_ic.H>
#include <mfix_bc.H>
#include <mfix_eb_if.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_reactions.H>
#include <mfix_leveldata.H>
#include <mfix_stepdata.H>
#include <mfix_mlmg_options.H>
#include <mfix_restart.H>
#include <mfix_run_on.H>
#include <mfix_rw.H>
#include <mfix_mass_balance.H>
#include <mfix_tag.H>
#include <mfix_timer.H>
#include <mfix_cmi.H>
#include <mfix_coupling.H>
#include <mfix_porous_media.H>
#include <mfix_eb.H>
#include <mfix_load_balance.H>
#include <mfix_avg_pc_parms.H>

// Mfix run type
enum struct RunType { Standard, Restart, PIC2DEM };

class mfix : public AmrCore
{
  public:
    static constexpr Real covered_val = 1.e40;

    // Constructor
    mfix ();

    // Destructor
    ~mfix ();


    // Initial Conditions
    MFIXInitialConditions m_initial_conditions;

    MFIXInitialConditions& ics () noexcept { return m_initial_conditions; }

    // Boundary conditions
    MFIXBoundaryConditions m_boundary_conditions;

    MFIXBoundaryConditions& bcs () noexcept { return m_boundary_conditions; }

    // Embedded boundaries parameters (flow, temperature, etc.)
    MFIXEmbeddedBoundaries m_embedded_boundaries;

    MFIXEmbeddedBoundaries& eb_parms () noexcept { return m_embedded_boundaries; }

    // Porous media regions
    MFIXPorousMedia m_porous_media;

    MFIXPorousMedia& pms () noexcept { return m_porous_media; }

    // Timer
    MFIXTimer m_timer;

    const MFIXTimer& timer () const { return m_timer; }
    MFIXTimer& timer () { return m_timer; }

    // Regions
    MFIXRegions regions;

    // Chemical species
    MFIXSpecies species;

    // Chemical reactions
    MFIXReactions reactions;

    // Fluid phase
    MFIXFluidPhase fluid;

    // Solids phase
    MFIXSolidsPhase solids;

    // DEM parameters
    MFIXDEM m_dem;

    // PIC parameters
    MFIXPIC m_pic;

    // The chemistry management interface class
    ChemistryManagementInterface m_cmi;

    // IO Class
    std::shared_ptr<MFIXMassBalance> m_mass_balance;

    // IO Class
    MFIXReadWrite* m_rw;

    // EB Class
    MFIXEB* m_eb;

    MFIXEB* eb () noexcept { return m_eb; }

    std::unique_ptr<LoadBalance> m_loadbalance;

    // Tags for AMR

    MFIXTag m_tag;

    MFIXTag& tag () noexcept { return m_tag; }

    // Level data
    MFIXLevelData& leveldata () noexcept { return m_level_data; }

    MFIXLevelData const& leveldata_const () const noexcept { return m_level_data; }

    LevelData* leveldata ( int const a_lev ) noexcept {
      AMREX_ASSERT( a_lev < m_level_data.m_level_data.size() );
      return m_level_data.m_level_data[a_lev].get();
    }

    // Particle container
    MFIXParticleContainer* pc;

    // Averaged particle data parameters
    MFIXAvgParticleParms m_avg_pc_parms;

    MFIXAvgParticleParms& avg_pc_parms () { return m_avg_pc_parms; }

    // MFIX class functions
    // ---------------------------------------------------------------------------//

    void InitParams ();

    void Init (amrex::Real time);

    BCList const& get_bc_list () const { return bc_list; }

    void InitParticlesFromScratch ();

    void init_advection ();

    void InitialRedistribution (amrex::Real l_time);

    void PostProjectionRedistribution (amrex::Real l_time, amrex::Real l_dt,
                                       const amrex::Vector<amrex::MultiFab*>& a_sigma);

    void PreProjectionRedistribution (amrex::Real l_time);

    void PostProjectionDiagnostics ( amrex::Real l_time,
                                     amrex::Vector< amrex::MultiFab      *> const& vel_g_in,
                                     amrex::Vector< amrex::MultiFab      *> const& p_g_in,
                                     amrex::Vector< amrex::MultiFab      *> const& gp_in,
                                     bool proj_for_small_dt);

    void Restart (std::string& restart_chkfile,
                  int& nstep,
                  amrex::Real& dt,
                  amrex::Real& time,
                  amrex::IntVect& Nrep);

    // Read (and replicate) grid data from checkpoint file
    void restart_leveldata_var ( amrex::IntVect a_Nrep,
                                 amrex::MultiFab* a_levMF,
                                 std::string const& a_prefix);

    // Create a temp particle container to read in checkpoint file then
    // copy to the runtime particle container.
    template <typename PatchType>
    void restart_from_old_chkpt ( const std::string& a_restart_file,
                                  const std::string& a_name,
                                  PatchType& a_patch );


    void average_thermodynamic_pressure ( const int lev,
                                          amrex::Real& thermo_p) const;

    void Regrid ( int const a_step );
    void Regrid ();

    void Evolve ();

    void EvolveSteadyState ();
    void Evolve_MOL ();
    void Evolve_Godunov ();

    void Evolve_OLD (); // To be removed

    void prepareToStep ( MFIXStepData& a_stepData );

    int m_use_new_depdt_algo = 0;

    bool has_particles () const noexcept
    { return (m_dem.solve() || m_pic.solve()); }


    // flag enabling level-set restart (i.e. prevent make_eb_* from rebuilding
    // the level-set data).
    bool levelset_restart = false;

    void compute_dt (amrex::Vector< amrex::MultiFab const* > const& a_avg_particle_data);

    void average_pc_data_to_grid ( amrex::Vector< MultiFab* > const& a_avg_data );

    void calc_txfr_particle ( amrex::Real const time, Real const dt,
                              amrex::Vector< amrex::MultiFab      *> const& a_vel,
                              amrex::Vector< amrex::MultiFab      *> const& a_T,
                              amrex::Vector< amrex::MultiFab const*> const& a_gradp,
                              amrex::Vector< amrex::MultiFab const*> const& a_divtau,
                              amrex::Vector< amrex::MultiFab const*> const& a_DufDt );

    void add_vel_src_implicit ( amrex::Real const a_dt,
                                amrex::Vector<amrex::MultiFab      *> const& a_vel,
                                amrex::Vector<amrex::MultiFab const*> const& a_epf,
                                amrex::Vector<amrex::MultiFab const*> const& a_rho,
                                amrex::Vector<amrex::MultiFab const*> const& a_S_p,
                                amrex::Vector<amrex::MultiFab const*> const& a_S_c);

    void add_enthalpy_txfr_implicit (amrex::Real dt,
                                     amrex::Vector<amrex::MultiFab*      > const& h_g_in,
                                     amrex::Vector<amrex::MultiFab*      > const& T_g_in,
                                     amrex::Vector<amrex::MultiFab const*> const& X_gk_in,
                                     amrex::Vector<amrex::MultiFab const*> const& txfr_in,
                                     amrex::Vector<amrex::MultiFab const*> const& rho_in,
                                     amrex::Vector<amrex::MultiFab const*> const& ep_g_in);

    void mfix_correct_small_cells (Vector<MultiFab*      > const& vel_in,
                                   Vector<MultiFab const*> const& ep_u_mac,
                                   Vector<MultiFab const*> const& ep_v_mac,
                                   Vector<MultiFab const*> const& ep_w_mac,
                                   Vector<MultiFab const*> const& eb_vel);

    static RunType m_run_type;

    bool run_type ( RunType a_run_type ) const
    { return (m_run_type == a_run_type); }

    void set_run_type ( RunType a_run_type ) noexcept
    { m_run_type = a_run_type; }

    CouplingOp m_coupling;

    amrex::Real deposit_volume_to_grid ( Real const a_time );
    void deposit_forces_to_grid ( Real const a_dt );

    void ComputeVariableFilter ( amrex::Vector< amrex::MultiFab const*> const& a_eps );

    ////////////////////////////////////////////////////////////////////////////
    //
    // Convection
    //
    ////////////////////////////////////////////////////////////////////////////

    void set_MAC_velocity_bcs (int a_lev,
                               amrex::MultiFab const* a_mac_rhs,
                               amrex::Array< amrex::MultiFab*,3 > a_mac_vel);

    void compute_MAC_projected_velocities (
             amrex::Real const a_time, amrex::Real const a_dt,
             amrex::Vector< amrex::MultiFab const*> const& a_vel,
             amrex::Vector< amrex::MultiFab const*> const& a_mac_rhs,
             amrex::Vector< amrex::Array<amrex::MultiFab      *,3>> const& a_mac_vel,
             amrex::Vector< amrex::MultiFab      *> const& a_mac_phi,
             amrex::Vector< amrex::MultiFab const*> const& a_epf,
             amrex::Vector< amrex::MultiFab const*> const& a_rho,
             amrex::Vector< amrex::MultiFab const*> const& a_txfr,
             amrex::Vector< amrex::MultiFab const*> const& a_eb_vel,
             amrex::Vector< amrex::MultiFab      *> const& a_vel_forces,
             amrex::Vector< amrex::MultiFab const*> const& a_divtau);

    ////////////////////////////////////////////////////////////////////////////
    //
    // Diffusion
    //
    ////////////////////////////////////////////////////////////////////////////

    void clear_solvers ();

    std::unique_ptr<MFIXDiffOpVelocity> m_diffOpVel;
    std::unique_ptr<MFIXDiffOpSpecies>  m_diffOpSpecies;
    std::unique_ptr<MFIXDiffOpEnergy>   m_diffOpEnergy;
    std::unique_ptr<MFIXDiffOpTracer>   m_diffOpTracer;
    std::unique_ptr<MFIXDiffOpScalar>   m_diffOpVoidFrac;
    std::unique_ptr<MFIXDiffOpScalar>   m_diffOpTxfr;

    MFIXDiffOpVelocity* diffOpVel () {

      if (!m_diffOpVel) {

        m_diffOpVel = std::make_unique<MFIXDiffOpVelocity>(
            nlev(), geom, grids, dmap, eb()->factory_const(),
            fluid, avg_pc_parms(), bcs(), bcs().get_velocity_bcrec() );

        if (eb_parms().has_flow()) { m_diffOpVel->setEBDirichlet(); }
      }
      return m_diffOpVel.get();
    }

    MFIXDiffOpSpecies* diffOpSpecies () {

      if (!m_diffOpSpecies) {
        m_diffOpSpecies = std::make_unique<MFIXDiffOpSpecies>(
            nlev(), geom, grids, dmap, eb()->factory_const(), fluid.nspecies(),
            fluid, bcs(), bcs().get_species_bcrec() );
      }
      return m_diffOpSpecies.get();
    }

    MFIXDiffOpEnergy* diffOpEnergy () {

      if (!m_diffOpEnergy) {
        m_diffOpEnergy = std::make_unique<MFIXDiffOpEnergy>(
            nlev(), geom, grids, dmap, eb()->factory_const(),
            fluid, bcs(), bcs().get_enthalpy_bcrec() );

        if (eb_parms().has_temperature()) { m_diffOpEnergy->setEBDirichlet(); }
      }
      return m_diffOpEnergy.get();
    }

    MFIXDiffOpTracer* diffOpTracer () {

      if (!m_diffOpTracer) {

        int const ncomp( fluid.ntracer() );

        m_diffOpTracer = std::make_unique<MFIXDiffOpTracer>(
            nlev(), geom, grids, dmap, eb()->factory_const(), ncomp, ncomp,
            fluid, bcs(), bcs().get_tracer_bcrec() );

        m_diffOpTracer->setDiffCoeff( fluid.tracer_diff_coeff() );
      }
      return m_diffOpTracer.get();
    }

    MFIXDiffOpScalar* diffOpVoidFrac () {

      if (!m_diffOpVoidFrac) {

        int const ncomp( 1 );

        m_diffOpVoidFrac = std::make_unique<MFIXDiffOpScalar>(
            nlev(), geom, grids, dmap, eb()->factory_const(), ncomp, ncomp,
            fluid, bcs(), bcs().get_volfrac_bcrec() );

        DepositionFilter const* depop_filter = m_coupling.getDepositionFilter();

        if ( depop_filter->ConstantSize() ) {
          Real const filter_size( depop_filter->size() );
          m_diffOpVoidFrac->setDiffCoeff(filter_size);
        }
      }
      return m_diffOpVoidFrac.get();
    }

    MFIXDiffOpScalar* diffOpTxfr () {

      if (!m_diffOpTxfr) {

        int const ncomp( leveldata().txfr(0)->nComp() );

        m_diffOpTxfr = std::make_unique<MFIXDiffOpScalar>(
            nlev(), geom, grids, dmap, eb()->factory_const(), ncomp, ncomp,
            fluid, bcs(), bcs().get_volfrac_bcrec() );

        DepositionFilter const* depop_filter = m_coupling.getDepositionFilter();

        if ( depop_filter->ConstantSize() ) {
          Real const filter_size( depop_filter->size() );
          m_diffOpTxfr->setDiffCoeff(filter_size);
        }
      }
      return m_diffOpTxfr.get();
    }

    ////////////////////////////////////////////////////////////////////////////
    //
    // Nodal Projection
    //
    ////////////////////////////////////////////////////////////////////////////

    std::unique_ptr<Hydro::NodalProjector> nodal_projector;

    ////////////////////////////////////////////////////////////////////////////
    //
    // PIC solids
    //
    ////////////////////////////////////////////////////////////////////////////

    void EvolveParcels (amrex::Real dt,
                        amrex::Real time,
                        amrex::RealVect& gravity,
                        const int ls_refinement,
                        LoadBalance* const a_loadbalance);

    void PostEvolveParcels ();

    void pic_iteration (const bool apply_forces,
                        const bool update_parcels,
                        const bool use_taylor_approx,
                        const Real advance_vel_p,
                        amrex::Real dt,
                        amrex::RealVect& gravity,
                        amrex::Vector< amrex::Array<amrex::MultiFab*,3> >& vel_s_in,
                        amrex::Vector< amrex::MultiFab* >& ep_s_out,
                        amrex::EBFArrayBoxFactory* ebfactory,
                        const int ls_refinement,
                        const MultiFab* ls_phi);

    void MFIX_CalcAvgSolidsVel (Vector< Array<MultiFab*,3> >& vel_s,
                                const bool do_deposition=false);

    void MFIX_CalcSolidsStress (amrex::Vector< amrex::MultiFab* >& ep_s_in,
                                amrex::RealVect& gravity_in,
                                LoadBalance* const a_loadbalance);

    void mfix_set_pic_velocity_bcs (Vector< MultiFab* > const& vel) const;

    void set_pic_velocity_bcs (const int lev,
                               FArrayBox& vel_fab,
                               const Box& domain) const;

    //! Tagging cells for refinement
    virtual void ErrorEst (int lev,
                           TagBoxArray & tags,
                           Real time,
                           int ngrow) override;

  protected:

    void set_gp0 ();

    static RealVect gp0;
    static RealVect gravity;

    BCList bc_list;

    void project_velocity ();

    void initial_iterations ();

  public:

    void apply_nodal_projection ( amrex::Real a_time, amrex::Real a_dt,
                                  amrex::Real a_prev_dt, bool proj_2,
                                  amrex::Vector< amrex::MultiFab      *> const& a_rhs,
                                  amrex::Vector< amrex::MultiFab      *> const& a_vel_old,
                                  amrex::Vector< amrex::MultiFab      *> const& a_vel_new,
                                  amrex::Vector< amrex::MultiFab      *> const& a_pert_p,
                                  amrex::Vector< amrex::MultiFab      *> const& a_grad_p,
                                  amrex::Vector< amrex::MultiFab const*> const& a_epf,
                                  amrex::Vector< amrex::MultiFab const*> const& a_rho,
                                  amrex::Vector< amrex::MultiFab const*> const& a_eb_vel);

    void compute_tra_forces (amrex::Vector<amrex::MultiFab      *> const& tra_forces,
                             amrex::Vector<amrex::MultiFab const*> const& density);

    void compute_vel_forces (amrex::Vector<amrex::MultiFab*      > const& vel_forces,
                             amrex::Vector<amrex::MultiFab const*> const& velocity,
                             amrex::Vector<amrex::MultiFab const*> const& density,
                             amrex::Vector<amrex::MultiFab const*> const& txfr_in,
                             bool include_pressure_gradient = true);

    void compute_vel_forces_on_level ( int lev,
                                       amrex::MultiFab& vel_forces,
                                       const amrex::MultiFab& velocity,
                                       const amrex::MultiFab& density,
                                       const amrex::MultiFab& txfr_in,
                                       bool include_pressure_gradient = true);

  protected:

    void avgDown (int crse_lev, const MultiFab & S_fine, MultiFab & S_crse);

  public:

    void MakeNewLevelFromScratch (int lev, Real time, const BoxArray & new_grids,
                                  const DistributionMapping & new_dmap) override;


    //! Apply MFIX-specific changes the to base grids.
    //!
    //! This function is only called by MakeNewGrids after computing a box array for
    //! the coarsest level and before calling MakeNewLevelFromScratch.
    virtual void PostProcessBaseGrids( amrex::BoxArray& a_grids) const override;


    //! Make a new level using provided BoxArray and DistributionMapping and
    //! fill with interpolated coarse level data. Called by AmrCore::regrid.
    void MakeNewLevelFromCoarse (int a_lev, amrex::Real a_time,
                                 BoxArray            const& a_grids,
                                 DistributionMapping const& a_dmap) override;

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse data. Called
    //! by AmrCore::regrid.
    void RemakeLevel ( int a_lev, amrex::Real a_time, const BoxArray & a_grids,
                       const DistributionMapping & a_dmap) override;

    //! Delete level data.  Called by AmrCore::regrid.
    void ClearLevel (int a_lev) override;

  private:

    void EvolveFluid ( MFIXStepData& a_step_data,
                       int nstep,
                       amrex::Real & dt,
                       amrex::Real & prev_dt,
                       const amrex::Real time,
                       amrex::Real stop_time,
                       amrex::Real & drag_timing);

    void EvolveSolids ( amrex::Real const a_dt,
                        amrex::Real const a_time,
                        int& a_nsubsetps,
                        amrex::Real& a_particle_timing,
                        amrex::Real& a_coupling_timing,
                        amrex::Vector< amrex::MultiFab* > const& a_DufDt,
                        amrex::Vector< amrex::MultiFab const*> const& a_avg_pc_data,
                        amrex::Vector< amrex::MultiFab      *> const& a_Teb);

    void fillpatch_mac (Vector< MultiFab* > const& ep_u_mac,
                        Vector< MultiFab* > const& ep_v_mac,
                        Vector< MultiFab* > const& ep_w_mac);

    int SteadyStateReached ( amrex::Real dt, int iter,
                             amrex::Vector< amrex::MultiFab const*> const& a_pg_prev);

  public:

    void compute_convective_term ( amrex::Vector< amrex::MultiFab      *> const& conv_u,
                                   amrex::Vector< amrex::MultiFab      *> const& conv_s,
                                   amrex::Vector< amrex::MultiFab      *> const& conv_X,
                                   amrex::Vector< amrex::MultiFab      *> const& vel_forces,
                                   amrex::Vector< amrex::MultiFab      *> const& tra_forces,
                                   amrex::Vector< amrex::MultiFab const*> const& vel_in,
                                   amrex::Vector< amrex::MultiFab const*> const& ep_g_in,
                                   amrex::Vector< amrex::MultiFab const*> const& ro_g_in,
                                   amrex::Vector< amrex::MultiFab const*> const& h_g_in,
                                   amrex::Vector< amrex::MultiFab const*> const& trac_in,
                                   amrex::Vector< amrex::MultiFab const*> const& X_gk_in,
                                   amrex::Vector< amrex::MultiFab const*> const& txfr_in,
                                   amrex::Vector< amrex::MultiFab const*> const& eb_vel,
                                   amrex::Vector< amrex::MultiFab const*> const& eb_scalars,
                                   amrex::Vector< amrex::MultiFab const*> const& eb_species,
                                   amrex::Vector< amrex::MultiFab      *> const& u_mac,
                                   amrex::Vector< amrex::MultiFab      *> const& v_mac,
                                   amrex::Vector< amrex::MultiFab      *> const& w_mac,
                                   Real l_dt, Real time);

    void apply_predictor ( MFIXStepData& a_stepData, amrex::Real time,
                           amrex::Real dt, amrex::Real prev_dt, bool proj_2);


    void apply_corrector ( MFIXStepData& a_step_Data, amrex::Real time,
                           amrex::Real dt, amrex::Real prev_dt, bool proj_2);

    void compute_MAC_proj_RHS ( amrex::Vector< amrex::MultiFab      * > const& a_rhs,
                                amrex::Vector< amrex::MultiFab const* > const& a_depdt,
                                amrex::Vector< amrex::MultiFab const* > const& a_epf,
                                amrex::Vector< amrex::MultiFab const* > const& a_rho,
                                amrex::Vector< amrex::MultiFab const* > const& a_Xf,
                                amrex::Vector< amrex::MultiFab const* > const& a_Tf,
                                amrex::Vector< amrex::MultiFab const* > const& a_lap_Tf,
                                amrex::Vector< amrex::MultiFab const* > const& a_div_hJ,
                                amrex::Vector< amrex::MultiFab const* > const& a_div_J,
                                amrex::Vector< amrex::MultiFab      * > const& a_txfr,
                                amrex::Vector< amrex::MultiFab const* > const& a_chem,
                                amrex::Real const& a_therm_p,
                                amrex::Real      & a_RHS_a_therm_p);

    void compute_nodal_proj_RHS ( amrex::Vector< amrex::MultiFab      * > const& a_rhs,
                                  amrex::Vector< amrex::MultiFab const* > const& a_depdt,
                                  amrex::Vector< amrex::MultiFab const* > const& a_epf_o,
                                  amrex::Vector< amrex::MultiFab const* > const& a_epf_n,
                                  amrex::Vector< amrex::MultiFab const* > const& a_rho_o,
                                  amrex::Vector< amrex::MultiFab const* > const& a_rho_n,
                                  amrex::Vector< amrex::MultiFab const* > const& a_hf_o,
                                  amrex::Vector< amrex::MultiFab const* > const& a_hf_n,
                                  amrex::Vector< amrex::MultiFab const* > const& a_Tf_n,
                                  amrex::Vector< amrex::MultiFab const* > const& a_dhdt,
                                  amrex::Vector< amrex::MultiFab const* > const& a_Xf_o,
                                  amrex::Vector< amrex::MultiFab const* > const& a_Xf_n,
                                  amrex::Vector< amrex::MultiFab const* > const& a_dXdt,
                                  amrex::Real const& a_therm_p,
                                  amrex::Real      & a_RHS_therm_p,
                                  amrex::Real const& a_dt);

    void idealgas_opensystem_rhs (Vector< MultiFab*       > const& rhs,
                                  Vector< MultiFab const* > const& enthalpy_rhs,
                                  Vector< MultiFab const* > const& species_rhs,
                                  Vector< MultiFab const* > const& ro_g,
                                  Vector< MultiFab const* > const& T_g,
                                  Vector< MultiFab const* > const& X_gk);

    void idealgas_closedsystem_rhs (amrex::Vector< amrex::MultiFab*       > const& a_rhs,
                                    amrex::Vector< amrex::MultiFab const* > const& a_epf,
                                    amrex::Vector< amrex::MultiFab const* > const& a_rho,
                                    amrex::Vector< amrex::MultiFab const* > const& a_Tf,
                                    amrex::Vector< amrex::MultiFab const* > const& a_Xf,
                                    amrex::Real const& a_therm_p,
                                    amrex::Real      & a_RHS_therm_p);

  private:

    void PruneBaseGrids(BoxArray &ba) const;

    void ChopGrids (const Box& domain, BoxArray& ba, int target_size) const;

    MFIXLevelData m_level_data;

  public:

    int nlev () const { return (finest_level+1); }

    // thermodynamic pressure
    amrex::Real therm_p;
    amrex::Real therm_po;

  private:

    MFIXParticleContainer* getParticleContainer() const { return pc; }

    amrex::Vector<int> const& get_velocity_iconserv () const noexcept { return m_iconserv_velocity; }
    amrex::Vector<int> const& get_density_iconserv () const noexcept { return m_iconserv_density; }
    amrex::Vector<int> const& get_enthalpy_iconserv () const noexcept { return m_iconserv_enthalpy; }
    amrex::Vector<int> const& get_tracer_iconserv () const noexcept { return m_iconserv_tracer; }
    amrex::Vector<int> const& get_species_iconserv () const noexcept { return m_iconserv_species; }
    //
    int const* get_velocity_iconserv_device_ptr () const noexcept {
      return m_iconserv_velocity_d.data(); }
    int const* get_density_iconserv_device_ptr () const noexcept {
      return m_iconserv_density_d.data(); }
    int const* get_enthalpy_iconserv_device_ptr () const noexcept {
      return m_iconserv_enthalpy_d.data(); }
    int const* get_tracer_iconserv_device_ptr () const noexcept {
      return m_iconserv_tracer_d.data(); }
    int const* get_species_iconserv_device_ptr () const noexcept {
      return m_iconserv_species_d.data(); }

    bool test_tracer_conservation = false;

    // If false then we update velocity using the convective   form, u dot grad u
    // If  true then we update velocity using the conservative form, del dot (u u)
    bool m_advect_momentum = false;

    AdvectionType m_advection_type = AdvectionType::Undefined;

    std::string m_redistribution_type = "StateRedist";

    bool m_grid_pruning = false;

    // Redistribute the velocity field before the nodal projection
    bool m_redistribute_before_nodal_proj = true;

    // Redistribute the velocity field after the nodal projection
    bool m_redistribute_nodal_proj = false;

    // Redistribute after the initial nodal projection
    bool m_redistribute_after_initial_nodal_proj = true;

    // Threshold volfrac for correcting small cell velocity in the predictor and corrector
    Real m_correction_small_volfrac = 1.e-4;

    // Tolerances and max iterations for Newton solver
    amrex::Real newton_reltol = 1.e-8;
    amrex::Real newton_abstol = 1.e-8;
    int newton_maxiter = 500;

    // If using Godunov, default to PLM (not PPM)
    bool m_godunov_ppm = false;

    // If using Godunov, default to including viscous/diffusive
    //    terms in the forcing terms
    bool m_godunov_include_diff_in_forcing = true;

    // If using Godunov, when to add the forcing terms in
    //    the construction of the "trans" velocities
    bool m_godunov_use_forces_in_trans = false;

    // If using Godunov, if true then don't include the pressure gradient
    //    in the forcing term passed into the godunov routine; instead
    //    use gradient of mac phi which contains the full pressure
    bool m_godunov_use_mac_phi = false;

    enum struct DiffusionType {
            Undefined, Explicit, Crank_Nicolson, Implicit
    };
    // MOL: Explicit predictor / Crank_Nicolson corrector
    // Godunov: Implicit
    DiffusionType m_predictor_diff_type = DiffusionType::Undefined;
    DiffusionType m_corrector_diff_type = DiffusionType::Undefined;

    amrex::Real enthalpy_source = 0;
    int update_mass             = 1;
    int update_momentum         = 1;
    int update_enthalpy         = 1;


    // Options to control mfix verbosity level
    int m_verbose = 0;


    // ***************************************************************
    // Nodal solve
    // ***************************************************************

    // Verbosity and MLMG parameters are now ParmParse with "nodal_proj" in the
    // inputs file
    // Examples: nodal_proj.verbose = 1
    //           nodal_proj.bottom_verbose = 1
    //           nodal_proj.maxiter
    //           nodal_proj.bottom_maxiter
    //           nodal_proj.bottom_rtol
    //           nodal_proj.bottom_atol
    //           nodal_proj.bottom_solver
    // More info at "AMReX-Hydro/Projections/hydro_NodalProjector.cpp"
    // Control MLMG behavior
    std::unique_ptr<MfixUtil::MLMGOptions> nodalproj_options;

    // ***************************************************************
    // MAC solve
    // ***************************************************************
    // Verbosity and MLMG parameters are now ParmParse with "mac_proj" in the
    // inputs file
    // Examples: mac_proj.verbose = 1
    //           mac_proj.bottom_verbose = 1
    //           mac_proj.maxiter
    //           mac_proj.bottom_maxiter
    //           mac_proj.bottom_rtol
    //           mac_proj.bottom_atol
    //           mac_proj.bottom_solver
    // More info at "AMReX-Hydro/Projections/hydro_MacProjector.cpp"
    std::unique_ptr<Hydro::MacProjector> macproj;
    // Control MLMG behavior
    std::unique_ptr<MfixUtil::MLMGOptions> macproj_options;
    // ***************************************************************

    int m_regrid_int = -1;

    // If true then print the name of the routine we are in
    bool ooo_debug = false;


    // Options to control time stepping
    Real m_cfl = 0.5;

    Real m_scale_init_dt = 1.0;

    int agg_grid_size = -1; // agglomeration for GMG

    int removeOutOfRange = 1; // Flag to remove out-of-range particles

    int sort_particle_int = -1; // # steps to sort particles

#ifdef AMREX_USE_GPU
    bool reduceGhostParticles = true;
#else
    bool reduceGhostParticles = false;
#endif


    static EBSupport m_eb_support_level;

    Periodicity p0_periodicity;

    Real sum_vol_orig;

    Real init_value = 9.87654321e32;


    // Here we set the number of ghost nodes of the field arrays.
    // Since this is  a class attribute, it will be propagated to
    // all the member functions.
    //
    // DO NOT USE  A LOCAL VARIABLE WITH THE SAME NAME
    //

    IntVect particle_sorting_bin = IntVect::TheUnitVector();

    bool do_initial_proj    = true;
    int  m_initial_iterations = 3;

    /***************************************************************************
    *                                                                          *
    * Used for fillpatching ...                                                *
    *                                                                          *
    ***************************************************************************/

  private:
    // Number of ghost cells for field arrays.
    int nghost_state () const { return 4; }

    // For Godunov, we need 1 ghost cell in addition to the Box we are filling
    // For MOL    , we need 0 ghost cells
    int nghost_force () const { return (advection_type() == AdvectionType::MOL) ? 0 : 1; }

    int nghost_mac () const {
      if (!m_eb->Factory(0).isAllRegular()) return (advection_type() == AdvectionType::Godunov) ? 4 : 3;
      return (advection_type() == AdvectionType::Godunov) ? 1 : 0;
    }

    AdvectionType advection_type () const noexcept { return m_advection_type;}

   // Undefined, Explicit, Crank_Nicolson, Implicit
    DiffusionType predictor_diff_type () const noexcept { return m_predictor_diff_type; }
    DiffusionType corrector_diff_type () const noexcept { return m_corrector_diff_type; }


    bool need_divtau () const {
      return ( (m_godunov_include_diff_in_forcing ||
                predictor_diff_type() != DiffusionType::Implicit) );
    }

    amrex::Vector<int> m_iconserv_velocity;
    amrex::Gpu::DeviceVector<int> m_iconserv_velocity_d;
    amrex::Vector<int> m_iconserv_density;
    amrex::Gpu::DeviceVector<int> m_iconserv_density_d;
    amrex::Vector<int> m_iconserv_enthalpy;
    amrex::Gpu::DeviceVector<int> m_iconserv_enthalpy_d;
    amrex::Vector<int> m_iconserv_tracer;
    amrex::Gpu::DeviceVector<int> m_iconserv_tracer_d;
    amrex::Vector<int> m_iconserv_species;
    amrex::Gpu::DeviceVector<int> m_iconserv_species_d;


};

#endif
