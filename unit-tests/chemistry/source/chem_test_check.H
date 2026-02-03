#ifndef CHEM_TEST_CHECK_H_
#define CHEM_TEST_CHECK_H_

#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_pc.H>

class chem_test_check {

  public:

    chem_test_check ( int const a_solve_species,
                      int const a_solve_enthalpy,
                      int const a_nspecies,
                      int const a_is_mixture,
                      int const a_write_particles,
                      amrex::Real const a_vol,
                      const ThermoPropertyData& a_fluid_props,
                      const ThermoPropertyData& a_solids_props );

    ~chem_test_check () {}

    void eulerian01 ( amrex::Real const a_time,
                      amrex::Vector< amrex::MultiFab* > const& a_field_data );

    void write ( amrex::Real const a_time,
                 amrex::Vector< amrex::MultiFab* > const& a_field_data,
                 MFIXParticleContainer* a_particles);

    amrex::Real L1norm ( ) {
      return (m_steps > 0)
        ? m_sum_abs_error / static_cast<amrex::Real>(m_steps)
        : std::numeric_limits<amrex::Real>::max();
    }

  private:

    int const m_solve_species;
    int const m_solve_enthalpy;
    int const m_nspecies;
    int const m_is_mixture;

    int const m_write_fluid;
    int const m_write_particles;

    amrex::Real const m_vol;

    const ThermoPropertyData& m_fluid_props;
    const ThermoPropertyData& m_solids_props;

    int m_verbose;

    int m_steps;
    amrex::Real m_sum_abs_error;
    amrex::Real m_max_abs_error;

  public:
    amrex::Real get_cell_value ( amrex::MultiFab const* const& a_MF,
                                 const amrex::IntVect& a_cell,
                                 int const a_comp);
};
#endif
