#include <AMReX_ParmParse.H>

#include <post_mfix.H>

using namespace amrex;

Vector< Vector<Real> >
post_mfix::
compute_averages (Vector<std::string> a_vars)
{

  int const ncomp = static_cast<int>(a_vars.size());

  Vector< MultiFab*> solids_ptr = m_solids->get_data();
  Vector< MultiFab*> fluid_ptr  = m_fluid->get_data();

  Vector< Vector<Real> > avgs;

  for (int lev(0); lev<get_nlev(); ++lev) {

    Real const grid_pts(static_cast<Real>((geom[lev].Domain()).numPts()));

    // Compute the global and phase averages for all vars

    MultiFab var_MF(grids[lev], dmap[lev], 1, 0);

    Vector<Real>   lev_avgs(ncomp, 0.0);

    // Compute global averages
    for( int comp(0); comp<ncomp; ++comp) {

      fill_var_MF(lev, a_vars[comp], &var_MF);

      lev_avgs[comp] = var_MF.sum(0, true);
    }

    ParallelDescriptor::ReduceRealSum(lev_avgs.data(), lev_avgs.size());

    for (int comp(0); comp < ncomp; ++comp)
    { lev_avgs[comp] /= grid_pts; }

    avgs.push_back(lev_avgs);

  }

  return avgs;

}


Vector< Vector<Real> >
post_mfix::
compute_fluctuations (Vector<std::string> a_vars,
                      Vector<Vector<Real>> a_avgs)
{
  int const ncomp = static_cast<int>(a_vars.size());

  Vector< MultiFab*> solids_ptr = m_solids->get_data();
  Vector< MultiFab*> fluid_ptr  = m_fluid->get_data();

  // Compute global fluctuations.

  Vector< Vector<Real> > flcts;

  for (int lev(0); lev<get_nlev(); ++lev) {

    Real const grid_pts(static_cast<Real>((geom[lev].Domain()).numPts()));

    MultiFab var_MF(grids[lev], dmap[lev], 1, 0);

    Vector<Real> lev_flcts(ncomp, 0.0);

    for( int comp(0); comp<ncomp; ++comp) {

      fill_var_MF(lev, a_vars[comp], &var_MF);

      // Compute global averages
      lev_flcts[comp] = calc_fluct(&var_MF, 0, a_avgs[lev][comp]);
    }

    ParallelDescriptor::ReduceRealSum( lev_flcts.data(), lev_flcts.size());

    for (int comp(0); comp < ncomp; ++comp) {
      lev_flcts[comp] = std::sqrt(lev_flcts[comp] / grid_pts);
    }

    flcts.push_back(lev_flcts);
  }
  return flcts;
}
