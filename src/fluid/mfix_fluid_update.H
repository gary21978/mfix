#ifndef MFIX_FLUID_UPDATE_H_
#define MFIX_FLUID_UPDATE_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <mfix_bc.H>

#include <mfix_diffusion_op.H>
#include <mfix_fluid.H>
#include <mfix_leveldata.H>
#include <mfix_stepdata.H>
#include <mfix_reactions.H>

enum struct StepType { Predictor, Corrector };

class FluidUpdate {

  public:

    FluidUpdate ( int const a_nlev, StepType const a_step_type,
                  amrex::Vector<amrex::Geometry> const& a_geom,
                  amrex::Vector<amrex::BoxArray> const& a_ba,
                  const MFIXFluidPhase& a_fluid,
                  MFIXLevelData& a_leveldata,
                  MFIXStepData& a_stepdata,
                  MFIXBoundaryConditions& a_bcs,
                  int const a_include_chem,
                  amrex::Real const a_time,
                  amrex::Real const a_dt);

    ~FluidUpdate () = default;

    void setExplicitDiffusion () { m_explicit_diffusion = 1; }

    int use_species_advection ( ) const noexcept
    { return (m_fluid.isMixture() && m_fluid.solve_species()); }

    void Density ( );

    void Species ( MFIXDiffOpSpecies* a_diffOpSpecies );

    void Energy ( amrex::Real const& a_pT_old, amrex::Real      & a_pT_new,
                  MFIXDiffOpEnergy* a_diffOp, amrex::Real const a_abstol,
                  amrex::Real const a_reltol, int const a_maxiter );

    void Tracers ( int const a_nscalars,
                   MFIXDiffOpTracer* a_diffOp);

    void Velocity ( int const a_advect_momentum,
                    int const a_include_src,
                    int const a_drag_includes_divtau,
                    MFIXDiffOpVelocity* a_diffOp);

    void VelocityImplicitUpdate ( amrex::MultiFab*       a_vel,
                                  amrex::MultiFab const* a_epf,
                                  amrex::MultiFab const* a_rho,
                                  amrex::MultiFab const* a_S_p,
                                  amrex::MultiFab const* a_S_c);

    void ComputeVelTxfrSrc ( AdvectionType const a_advection_type,
                             int const a_advect_momentum,
                             int const a_include_vm);

    void ComputeAcceleration ( AdvectionType const a_adv_type,
                               int const a_advect_momentum);

  protected:

    int nlev () { return m_nlev; }

    int finestLevel () { return m_nlev-1; }

    StepData& predictor () const noexcept { return m_stepdata.predictor(); }
    StepData& corrector () const noexcept { return m_stepdata.corrector(); }

    MFIXLevelData& leveldata() const noexcept { return m_leveldata; }
    MFIXStepData& stepdata() const noexcept { return m_stepdata; }

  private:

    int const m_nlev;

    StepType const m_step_type;

    amrex::Vector<amrex::Geometry> const& m_geom;

    amrex::Vector<amrex::BoxArray> const& m_ba;

    const MFIXFluidPhase& m_fluid;

    MFIXLevelData& m_leveldata;

    MFIXStepData& m_stepdata;

    MFIXBoundaryConditions& m_bcs;

    int const m_include_chem;

    amrex::Real const m_time;

    amrex::Real const m_dt;

    InterphaseChemTxfrIndexes m_chem_idxs;

    int m_explicit_diffusion = 0;

    void average_down ( int const a_lev, amrex::Vector<amrex::MultiFab*> a_MF) {

      amrex::IntVect ref_ratio = m_geom[a_lev+1].Domain().size() /
                                 m_geom[a_lev  ].Domain().size();

      amrex::EB_average_down( *a_MF[a_lev+1], *a_MF[a_lev], 0,
          a_MF[a_lev]->nComp(), ref_ratio);
    }

};
#endif
