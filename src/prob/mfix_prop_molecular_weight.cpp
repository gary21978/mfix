#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>

#include <mfix_prop_molecular_weight.H>
#include <mfix_reporter.H>

using namespace amrex;

int MolecularWeight::
define  ( std::string base )
{
  ParmParse ppMFIX("mfix");
  ParmParse ppSpecies("species");
  ParmParse ppBase(base); // e.g., fluid0

  // Check if solving species and if so, get the list defining base
  int solve_species(0);
  ppMFIX.query("solve_species", solve_species);

  Vector<std::string> species;
  if (solve_species) {
    ppBase.queryarr("species", species);
  }

  m_ncomp = amrex::max(1, static_cast<int>(species.size()));

  m_h_MW.resize(m_ncomp, 0.);
  m_d_MW.resize(m_ncomp, 0.);

  int readMW(0);

  if (solve_species) {

    for( int n(0); n<m_ncomp; ++n) {
      readMW += ReadEntry("species."+species[n], n);
    }
  } else {

    readMW += ReadEntry(base, 0);

  }

  if (readMW == m_ncomp) {

    m_defined = 1;

    // All checks passed, copy from host to device.
    Gpu::copyAsync(Gpu::hostToDevice, m_h_MW.begin(),
        m_h_MW.end(), m_d_MW.begin());
  }

  return 0;
}

int MolecularWeight::
ReadEntry ( std::string base, int const comp)
{
  ParmParse pp(base);

  return (pp.query("molecular_weight", m_h_MW[comp]));
}
