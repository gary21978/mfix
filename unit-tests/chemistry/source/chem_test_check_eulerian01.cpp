#include <chem_test_check.H>

using namespace amrex;

void
chem_test_check::
eulerian01 (Real const a_time,
            Vector< MultiFab* > const& a_field_data )
{
#ifndef CHEM_TEST_EULERIAN01
  ignore_unused(a_time);
  ignore_unused(a_field_data);
  m_steps = -1;
  m_max_abs_error = std::numeric_limits<amrex::Real>::max();
  m_sum_abs_error = std::numeric_limits<amrex::Real>::max();
#else

  DualGridAuxIndexes fluid_idxs(m_solve_enthalpy, m_nspecies, m_is_mixture);

  int const lev(0);

  Real abserr = std::abs( ( 0.5 - 0.25*std::sin(6.*M_PI*a_time)) -
    get_cell_value(a_field_data[lev], {0,0,0}, fluid_idxs.X_gk));

  m_steps += 1;
  m_sum_abs_error += abserr;
  m_max_abs_error = amrex::max(m_max_abs_error, abserr);
#endif
}
