#include <AMReX.H>
#include <AMReX_Print.H>

#include <fml_reporter.H>
#include <FilterML.H>

using namespace amrex;

// Destructor
FilterML::
~FilterML ()
{
  for (auto& lev_alpha_p : m_alpha_p) { lev_alpha_p.reset( nullptr ); }
  for (auto& lev_alpha_f : m_alpha_f) { lev_alpha_f.reset( nullptr ); }

  if ( m_fluid     ) { delete m_fluid;     }
  if ( m_particles ) { delete m_particles; }
  if ( m_solids    ) { delete m_solids;    }
}

// Constructor
FilterML::
FilterML ( std::string a_plotfile )
  : m_verbose(0)
  , m_plotfile(a_plotfile)
{
  fml_plotfile pf(a_plotfile);

  m_max_level = pf.finestLevel();

  int const nlev = get_nlev();

  m_geom.resize(nlev);
  m_grids.resize(nlev);
  m_dmap.resize(nlev);
  m_ref_ratio.resize(nlev);

  RealBox rb(pf.probLo(), pf.probHi());

  for (int lev(0); lev<nlev; ++lev) {

    Box dom = pf.probDomain(lev);

    m_grids[lev] = pf.boxArray(lev);
    m_dmap[lev] = pf.DistributionMap(lev);

    m_geom[lev].define(dom, rb, 0, {1,1,1});

    m_ref_ratio[lev] = pf.refRatio(lev);
  }

  { // Read and store fluid data
    Print() << "\nLoading grid data\n";
    m_fluid = new fml_fluid(m_max_level, m_geom, m_dmap, m_grids, &pf);
  }

  if (pf.has_particles() ) {
    Print() << "Loading particle data\n";
    m_particles = new fml_particles(m_max_level, m_geom, m_dmap, m_grids, m_ref_ratio, &pf);
  }

}


void
FilterML::
create_Eulerian_solids (Vector<std::string> a_vars)
{

  if ( !m_particles->contains("volume") && !m_particles->contains("radius") ) {
    fml_report::Error(__FILE__, __LINE__) << " Error:"
      << " Plotfile does not contain particle radius or volume.\n"
      << " Unable to proceed.";
  }

  // Convert particle data to continuum variables
  Print() << "\nCreating Eulerain solids\n";
  m_solids = new fml_solids(m_max_level, m_geom, m_dmap, m_grids, m_particles, a_vars);
  m_solids->reset();

  m_alpha_p.resize(get_nlev());
  m_alpha_f.resize(get_nlev());

  Print() << "\nComputing alpha_p and alpha_f\n";
  for (int lev(0); lev < get_nlev(); lev++) {

    m_alpha_p[lev].reset(new MultiFab(m_grids[lev], m_dmap[lev], 1, 1));
    m_alpha_p[lev]->setVal(0.0);
    { int const srccomp = m_solids->index("alpha");
      int const dstcomp = 0;
      int const numcomp = 1;
      int const nghost  = 0;

      MultiFab const* const src_ptr = m_solids->get_const_data(lev);
      MultiFab::Copy(*m_alpha_p[lev], *src_ptr, srccomp, dstcomp, numcomp, nghost);
    }

    m_alpha_f[lev].reset(new MultiFab(m_grids[lev], m_dmap[lev], 1, 1));
    m_alpha_f[lev]->setVal(1.0);
    { int const srccomp = 0;
      int const dstcomp = 0;
      int const numcomp = 1;
      int const nghost  = 0;

      MultiFab::Subtract(*m_alpha_f[lev], *m_alpha_p[lev], srccomp, dstcomp, numcomp, nghost);
    }
    m_alpha_p[lev]->FillBoundary(m_geom[lev].periodicity());
    m_alpha_f[lev]->FillBoundary(m_geom[lev].periodicity());
  }

}


void
FilterML::
fill_var_MF ( int const a_lev, std::string const a_var,
              MultiFab* a_MF, int const a_dstcomp)
{
  int const verbose(0);

  int const numcomp = 1;
  int const  nghost = 0;

  Vector<std::string> var_comps = get_var_comps(a_var);
  int const var_ncomp = static_cast<int>(var_comps.size());

  for (int n(0); n<var_ncomp; ++n) {

    std::string entry = var_comps[n];

    int srccomp;

    MultiFab const* src_ptr;

    a_MF->setVal(1.0, a_dstcomp, 1, a_MF->nGrow());

    if (entry.compare("alpha_p") == 0) {

      if (verbose) { Print() << "Found " << entry << "\n"; }

      srccomp = 0;
      src_ptr = get_const_alpha_p(a_lev);

    } else if (entry.compare("alpha_f") == 0) {

      if (verbose) { Print() << "Found " << entry << "\n"; }

      srccomp = 0;
      src_ptr = get_const_alpha_f(a_lev);

    } else if (m_solids->contains(entry)) {

      srccomp = m_solids->index(entry);
      src_ptr = m_solids->get_const_data(a_lev);

      if (verbose) { Print() << "Found " << entry << " in solids data.\n"; }

    } else if (m_fluid->contains(entry)) {

      srccomp = m_fluid->index(entry);
      src_ptr = m_fluid->get_const_data(a_lev);

      if (verbose) { Print() << "Found " << entry << " in fluid data.\n"; }

    } else {
      amrex::Print() << "Unable to locate " << entry << "!\n";
      amrex::Abort();
    }

    MultiFab::Multiply(*a_MF, *src_ptr, srccomp, a_dstcomp, numcomp, nghost);

  }
}
