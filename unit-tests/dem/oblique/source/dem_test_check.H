#ifndef DEM_TEST_CHECK_H_
#define DEM_TEST_CHECK_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_Vector.H>

#include <mfix_pc.H>
#include <mfix_dem.H>

class dem_test_check
{
  public:
    dem_test_check ();

    bool compare (MFIXParticleContainer& pc, int nlev) const;

    bool rolling_friction_inputs (MFIXDEM& dem) const;

  private:
    amrex::Real m_epsilon;
    amrex::Vector<amrex::RealVect> m_vel;
    amrex::Vector<amrex::RealVect> m_omega;
};

#endif
