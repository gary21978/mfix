#include <AMReX.H>
#include <AMReX_Print.H>

#include <post_mfix.H>

using namespace amrex;

post_mfix::
~post_mfix ()
{

  for (auto& lev_alpha_f : m_alpha_f) { lev_alpha_f.reset( nullptr ); }
  for (auto& lev_alpha_p : m_alpha_p) { lev_alpha_p.reset( nullptr ); }

  if ( m_fluid     ) { delete m_fluid; }
  if ( m_particles ) { delete m_particles; }
  if ( m_solids    ) { delete m_solids; }

}

// Constructor
post_mfix::
post_mfix ( amrex::Geometry const& a_level_0_geom,
            amrex::AmrInfo const& a_amr_info,
            pm_plotfile* a_plotfile)
  : amrex::AmrCore(a_level_0_geom, a_amr_info)
  , m_verbose(0)
{

  int const nlev = a_plotfile->get_nlev();
  Real const time = a_plotfile->get_time();

  for (int lev(0); lev < nlev; lev++) {

    const BoxArray& ba = a_plotfile->get_ba(lev);
    DistributionMapping dm = a_plotfile->get_dm(lev);

    MakeNewLevelFromScratch(0, time, ba, dm);
  }

  { // Read and store fluid data

    Print() << "\nLoading grid data\n";
    m_fluid = new pm_fluid(maxLevel(), Geom(), DistributionMap(),
                    boxArray(), a_plotfile );
  }


  { // Create particle container

    Print() << "Loading particle data\n";
    m_particles = new pm_particles(maxLevel(), Geom(), DistributionMap(),
                    boxArray(), a_plotfile);
    Print() << "Total number of particles: " << a_plotfile->particle_count() << "\n";
  }

  // Convert particle data to continuum variables
  Print() << "\nCreating Eulerain solids\n";
  m_solids = new pm_solids(maxLevel(), Geom(), DistributionMap(),
                    boxArray(), m_particles);
  m_solids->reset();


  m_alpha_p.resize(get_nlev());
  m_alpha_f.resize(get_nlev());

  Print() << "\nComputing alpha_p and alpha_f\n";
  for (int lev(0); lev < nlev; lev++) {

    m_alpha_p[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 1));
    m_alpha_p[lev]->setVal(0.0);
    { int const srccomp = m_solids->index("alpha");
      int const dstcomp = 0;
      int const numcomp = 1;
      int const nghost  = 0;

      MultiFab const* const src_ptr = m_solids->get_const_data(lev);
      MultiFab::Copy(*m_alpha_p[lev], *src_ptr, srccomp, dstcomp, numcomp, nghost);
    }

    m_alpha_f[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 1));
    m_alpha_f[lev]->setVal(1.0);
    { int const srccomp = 0;
      int const dstcomp = 0;
      int const numcomp = 1;
      int const nghost  = 0;

      MultiFab::Subtract(*m_alpha_f[lev], *m_alpha_p[lev], srccomp, dstcomp, numcomp, nghost);
    }
    m_alpha_p[lev]->FillBoundary(geom[lev].periodicity());
    m_alpha_f[lev]->FillBoundary(geom[lev].periodicity());
  }

}

void
post_mfix::
MakeNewLevelFromScratch (int a_lev, Real /*time*/,
                         const BoxArray& a_new_grids,
                         const DistributionMapping& a_new_dmap)
{
    if (m_verbose) {
      amrex::Print() << "MakeNewLevelFromScratch" << std::endl;

      if (a_lev == 0) {
        std::cout << "Making level 0 with box array\n" << a_new_grids << "\n";
      } else {
        Print() << "Setting refined region at level " << a_lev
                << " to " << a_new_grids << "\n";
      }
    }

    if (m_verbose > 1) {
    }

    SetBoxArray(a_lev, a_new_grids);
    SetDistributionMap(a_lev, a_new_dmap);

    if (m_verbose > 1) {
      amrex::Print() << "SETTING NEW GRIDS IN MAKE NEW LEVEL " << a_new_grids << std::endl;
      amrex::Print() << "SETTING NEW DMAP IN MAKE NEW LEVEL " << a_new_dmap << std::endl;
    }
}

void
post_mfix::
fill_var_MF ( int const a_lev, std::string const a_var,
              MultiFab* a_MF, int const a_dstcomp)
{
  a_MF->setVal(1.0, a_dstcomp, 1, a_MF->nGrow());

  int const verbose(0);

  //Print() << "\nFilling variable " << a_var << "\n";

  Vector<std::string> var_comps = get_var_comps(a_var);
  int const var_ncomp = static_cast<int>(var_comps.size());

  int const numcomp = 1;
  int const  nghost = 0;

  int srccomp;
  MultiFab const* src_ptr;

  for (int n(0); n<var_ncomp; ++n) {

    std::string entry = var_comps[n];

    if (is_gradient(entry)) {

      size_t first = entry.find("(");
      size_t last  = entry.rfind(")");

      std::string grad_var = entry.substr(first+1,last-first-1);

      MultiFab tmp_MF(grids[a_lev], dmap[a_lev], 1, 1);

      if (verbose) { Print() << "Computing gradient of " << grad_var << "\n";}
      int const dir = grad_dir(entry);
      compute_grad( a_lev, dir, grad_var, &tmp_MF);

      srccomp = 0;
      MultiFab::Multiply(*a_MF, tmp_MF, srccomp, a_dstcomp, numcomp, nghost);

    } else if (entry.find("sigma_") != std::string::npos) {

      if (verbose) { Print() << "Computing fluid stress component " << entry << "\n"; }

      MultiFab tmp_MF(grids[a_lev], dmap[a_lev], 1, 1);

      compute_stress( a_lev, entry, &tmp_MF);

      srccomp = 0;
      MultiFab::Multiply(*a_MF, tmp_MF, srccomp, a_dstcomp, numcomp, nghost);

    } else if (entry.find("sigPp_") != std::string::npos) {

      if (verbose) { Print() << "Computing granular viscosity term " << entry << "\n"; }

      MultiFab tmp_MF(grids[a_lev], dmap[a_lev], 1, 1);

      compute_granular_viscosity( a_lev, entry, &tmp_MF);

      srccomp = 0;
      MultiFab::Multiply(*a_MF, tmp_MF, srccomp, a_dstcomp, numcomp, nghost);

    } else {

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
}
