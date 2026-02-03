#ifndef MFIX_MASS_BALANCE_H
#define MFIX_MASS_BALANCE_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_leveldata.H>
#include <mfix_pc.H>
#include <mfix_bc.H>
#include <mfix_dem.H>
#include <mfix_pic.H>


class MFIXMassBalance
{
  public:
    MFIXMassBalance (amrex::Vector<amrex::Geometry> const& a_geom,
                     amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& a_ebfactory,
                     MFIXLevelData& a_leveldata,
                     MFIXParticleContainer* a_pc,
                     MFIXFluidPhase const& a_fluid,
                     MFIXSolidsPhase const& a_solids,
                     MFIXDEM const& a_dem,
                     MFIXPIC const& a_pic,
                     BCList const& a_bclist);

    void Initialize ();

    void set_pc (MFIXParticleContainer* a_pc) {
      m_pc = a_pc;
    }

    void ComputeMassAccum (const int report_mass_balance,
                           const int offset = 1);

    void ComputeMassProduction (amrex::Gpu::DeviceVector<amrex::Real> const& d_fluid_prod,
                                amrex::Gpu::DeviceVector<amrex::Real> const& d_solids_prod);

    void ComputeMassFlux (amrex::Vector< amrex::MultiFab const*> const& flux_x,
                          amrex::Vector< amrex::MultiFab const*> const& flux_y,
                          amrex::Vector< amrex::MultiFab const*> const& flux_z,
                          const int scomp,
                          const int ncomp,
                          const bool fluxes_are_area_weighted,
                          const int eb_has_flow,
                          amrex::Vector< amrex::MultiFab const*> const& eb_vel_in,
                          amrex::Vector< amrex::MultiFab const*> const& eb_species_in,
                          const amrex::Real dt);

    void ComputeMassOutflow (int const a_lev);

    void InitMassBalance (const int report_mass_balance,
                          const int nspecies);

    void ResetMassBalance (const int a_n);

    amrex::Vector<amrex::Real>& fluid_mass_accum () { return m_fluid_mass_accum; }
    amrex::Vector<amrex::Real>& fluid_mass_inflow () { return m_fluid_mass_inflow; }
    amrex::Vector<amrex::Real>& fluid_mass_outflow () { return m_fluid_mass_outflow; }
    amrex::Gpu::HostVector<amrex::Real>& fluid_mass_prod () { return m_fluid_mass_prod; }

    amrex::Vector<amrex::Real>& solids_mass_accum () { return m_solids_mass_accum; }
    amrex::Vector<amrex::Real>& solids_mass_inflow () { return m_solids_mass_inflow; }
    amrex::Vector<amrex::Real>& solids_mass_outflow () { return m_solids_mass_outflow; }
    amrex::Gpu::HostVector<amrex::Real>& solids_mass_prod () { return m_solids_mass_prod; }

    MFIXLevelData& leveldata () noexcept { return m_leveldata; }

  private:
    const int m_nlev;

    amrex::Vector<amrex::Real> m_fluid_mass_accum;
    amrex::Vector<amrex::Real> m_fluid_mass_inflow;
    amrex::Vector<amrex::Real> m_fluid_mass_outflow;
    amrex::Gpu::HostVector<amrex::Real> m_fluid_mass_prod;

    amrex::Vector<amrex::Real> m_solids_mass_accum;
    amrex::Vector<amrex::Real> m_solids_mass_inflow;
    amrex::Vector<amrex::Real> m_solids_mass_outflow;
    amrex::Gpu::HostVector<amrex::Real> m_solids_mass_prod;

//    amrex::Long m_total_numparticle;

    amrex::Vector<amrex::Geometry> const& m_geom;
    amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& m_ebfactory;
    MFIXLevelData& m_leveldata;
    MFIXParticleContainer* m_pc;
    MFIXFluidPhase const& m_fluid;
    MFIXSolidsPhase const& m_solids;
    MFIXDEM const& m_dem;
    MFIXPIC const& m_pic;
    BCList const& m_bclist;
};

#endif
