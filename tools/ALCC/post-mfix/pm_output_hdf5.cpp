#include <AMReX_ParmParse.H>

#include <post_mfix.H>

using namespace amrex;

void
post_mfix::
output_hdf5 (std::string a_plotfile)
{
#ifdef AMREX_USE_HDF5

  ParmParse pp_HDF5("hdf5");

  size_t last = a_plotfile.find_last_of("/\\");
  if (last == std::string::npos) { last = 0; }

  // Get the index of the current plot file.
  int const first = a_plotfile.find_first_of("0123456789",last);
  int nstep = std::stoi( a_plotfile.substr(first) );

  std::string path;
  if (!pp_HDF5.query("path",path) )
  { path = (last == 0) ? "./" : a_plotfile.substr(0,last+1); }

  if ( path.back() != '/' ) { path += "/"; }

  std::string hdf5_compression = "None@0";
#if (defined AMREX_USE_HDF5_ZFP) || (defined AMREX_USE_HDF5_SZ)
  pp_HDF5.query("compression", hdf5_compression);
#endif

  const std::string& hdf5_dir = amrex::Concatenate(path+"hdf5_",nstep);

#if 0
  Print() << "plotfile: " << a_plotfile << "\n"
          << "       path: " << path  << "\n"
          << "      nstep: " << nstep << "\n"
          << "       hdf5: " << hdf5_dir << "\n"
          << "compression: " << hdf5_compression << "\n";
#endif

  Vector<int> istep;
  istep.resize(get_nlev(),nstep);
  Real time = 0.;

  Vector<std::unique_ptr<MultiFab>> HDF5_MF(get_nlev());
  Vector<std::string> varnames;

  for (int lev(0); lev<get_nlev(); ++lev) {

    // Fluid level data
    auto const lev_data = m_fluid->get_const_data(lev);

    // Number of components (all fluid plus fluid volume fraction.
    int const ncomp = lev_data->nComp()+1;

    HDF5_MF[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, 0));

    MultiFab::Copy(*HDF5_MF[lev].get(), *get_const_alpha_f(lev), 0, 0, 1, 0);
    varnames.push_back("alpha_f");

    for (int comp(0); comp<lev_data->nComp(); ++comp) {
      MultiFab::Copy(*HDF5_MF[lev].get(), *lev_data, comp, comp+1, 1, 0);
      varnames.push_back(m_fluid->varName(comp));
    }
  }

  // Not writing particle real data?!
  Vector<int> write_real_comp(m_particles->nComps_real(),1);
  Vector<int> write_int_comp(m_particles->nComps_int(),1);

  Vector<std::string> real_comp_names = m_particles->VarNames_real();
  Vector<std::string> int_comp_names = m_particles->VarNames_int();

  AMREX_ALWAYS_ASSERT( write_real_comp.size() == real_comp_names.size());
  AMREX_ALWAYS_ASSERT( write_int_comp.size() == int_comp_names.size());

  // Create a directory to group grid and particle data.
  if (ParallelDescriptor::IOProcessor()) {
    if ( ! amrex::UtilCreateDirectory(hdf5_dir, 0755) )
    { amrex::CreateDirectoryFailed(hdf5_dir); }
  }
  ParallelDescriptor::Barrier();

  const std::string& hdf5_grid = hdf5_dir + "/grid_data";

  WriteMultiLevelPlotfileHDF5SingleDset(hdf5_grid, get_nlev(),
                                        GetVecOfConstPtrs(HDF5_MF), varnames,
                                        Geom(), time, istep, refRatio(),
                                        hdf5_compression);

  m_particles->WritePlotFileHDF5(hdf5_dir, "particles",
                                 write_real_comp, write_int_comp,
                                 real_comp_names, int_comp_names,
                                 hdf5_compression);

  ParallelDescriptor::Barrier();

  for (int lev(0); lev<get_nlev(); ++lev)
  { HDF5_MF[lev].reset( nullptr ); }

#endif
}
