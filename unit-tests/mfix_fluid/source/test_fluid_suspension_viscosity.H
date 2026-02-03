#ifndef TEST_FLUID_SUSPENSION_VISCOSITY_H_
#define TEST_FLUID_SUSPENSION_VISCOSITY_H_

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>

#include <memory>

class test_fluid_suspension_viscosity
{
  public:

    test_fluid_suspension_viscosity ();
    ~test_fluid_suspension_viscosity ();

    int einstein ();
    int brinkman ();
    int roscoe ();
    int cheng_law ();
    int sato ();

  private:

    void init_fabs (const amrex::Box& bx, amrex::Array4<amrex::Real>& epg,
                    amrex::Array4<amrex::Real>& rho);

    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::MultiFab* vel_fab;
    amrex::MultiFab* epf_fab;
    amrex::MultiFab* rho_fab;
    amrex::MultiFab* Tf_fab;
    amrex::MultiFab* avg_part_fab;
    std::unique_ptr<amrex::EBFArrayBoxFactory> ebfact;

    int const PASS = 0;
    int const FAIL = 1;
    amrex::Real const MOL_VISC_CONST = 1.2e-5;
    amrex::Real const DUMMY = 1.e30;
    amrex::Real BRINKMAN_CONST = 2.5;
    amrex::Real ROSCOE_C1 = 0.65;
    amrex::Real ROSCOE_C2 = 2.5;
    amrex::Real CHENG_LAW_CONST = 2.5;
    amrex::Real MAX_EFF_VISC_FAC = 1.e5;
    amrex::Real SATO_CONST = 0.7;

};
#endif
