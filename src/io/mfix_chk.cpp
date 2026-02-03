#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_AmrCore.H>
#include <AMReX_buildInfo.H>
#include <AMReX_Geometry.H>

#include <mfix_rw.H>
#include <mfix_fluid.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_versions.H>

using namespace amrex;

void MFIXReadWrite::
WriteCheckHeader ( const std::string& a_name, int a_nstep,
                   Real a_dt, Real a_time) const
{
  if (ParallelDescriptor::IOProcessor()) {

    std::string HeaderFileName(a_name + "/Header");
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    HeaderFile.open(HeaderFileName, std::ofstream::out   |
                    std::ofstream::trunc |
                    std::ofstream::binary);

    if ( ! HeaderFile.good() )
        amrex::FileOpenFailed(HeaderFileName);

    HeaderFile.precision(17);

    version::checkpoint version;

    HeaderFile << "Checkpoint version: " << version.current_str() << "\n";

    const int nlevels = nlev;
    HeaderFile << nlevels << "\n";

    // Time stepping controls
    HeaderFile << a_nstep << "\n";
    HeaderFile << a_dt << "\n";
    HeaderFile << a_time << "\n";

    // Geometry
    HeaderFile << RealVect(geom[0].ProbLo()) << "\n";
    HeaderFile << RealVect(geom[0].ProbHi()) << "\n";
    HeaderFile << geom[0].Domain().size() << "\n";

    Real small_volfrac(0.);
    ParmParse pp("eb2");
    pp.query("small_volfrac", small_volfrac);
    HeaderFile << small_volfrac << "\n";

    // BoxArray
    for (int lev = 0; lev < nlevels; ++lev)
    {
        grids[lev].writeOn(HeaderFile);
        HeaderFile << "\n";
    }
  }
}


void
MFIXReadWrite::
WriteCheckPointFile ( std::string& a_check_file, int a_nstep,
                      Real a_dt, Real a_time)
{
  BL_PROFILE("mfix::WriteCheckPointFile()");
  const std::string level_prefix {"Level_"};

  const std::string& checkpointname = amrex::Concatenate(
      a_check_file, a_nstep, m_min_digits );

  Print() << "\n\t Writing checkpoint " << checkpointname << '\n';

  const int nlevels = nlev;
  PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

  WriteCheckHeader(checkpointname, a_nstep, a_dt, a_time);

  WriteJobInfo(checkpointname);

  if (fluid.solve()) {

    for (int lev(0); lev < nlevels; ++lev) {

      // This writes all three velocity components
      VisMF::Write( *(leveldata().vel(lev)),
          MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "u_g"));

      // This writes all three pressure gradient components
      VisMF::Write( *(leveldata().grad_p(lev)),
          MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "gpx"));

      // volume fraction
      VisMF::Write( *(leveldata().epf(lev)),
          MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "ep_g"));

      // pressure
      VisMF::Write( *(leveldata().pert_p(lev)),
          MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "p_g"));

      // density
      VisMF::Write( *(leveldata().rho(lev)),
          MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "ro_g"));


      // temperature
      if (leveldata(lev)->has_temperature()) {
        VisMF::Write( *leveldata().T(lev),
            amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "T_g"));
      }

      // species mass fraction
      if (fluid.solve_species()) {
        VisMF::Write( *(leveldata().X(lev)),
            MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "X_gk"));
      }

      // thermodynamic pressure
      if (fluid.solve_enthalpy() && fluid.constraint.isIdealGasClosedSystem() ) {

        if (ParallelDescriptor::IOProcessor()) {

          std::string filename =
              MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "thermodynamic_p_g");

          std::ofstream filestream;
          filestream.open(filename, std::ios::out | std::ios::trunc);
          if (!filestream.good()) { amrex::FileOpenFailed(filename); }

          filestream << m_therm_p << "\n";
          filestream.close();
        }
      }

      // level set
      if ( m_dem.solve() || m_pic.solve() ) {
        VisMF::Write( *(level_sets[lev].get()),
            MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "level_sets"));
      }

    }
  }

  if ( m_dem.solve() || m_pic.solve() ) {
     pc->Checkpoint(checkpointname, "particles");
  }


  if (m_dem.solve() || m_pic.solve()) {

    // The level set might have a higher refinement than the mfix level.
    //      => Current mechanism for saving checkpoint files requires the
    //         same BoxArray for all MultiFabss on the same level
    // NOTE: the unrefined level-set (and the multi-level level set) are
    // both saved with the standard checkpoint file.
    std::stringstream raw_ls_name;
    raw_ls_name << checkpointname << "/ls_raw";

    // There is always a level 1 in the level_sets array:
    //    level_sets.size() == amrex::max(2, maxLevel())
    VisMF::Write( * level_sets[1], raw_ls_name.str() );

    // Also save the parameters necessary to re-build the LSFactory
    int levelset_params[] = { levelset_refinement,
                              levelset_pad,
                              levelset_eb_refinement,
                              levelset_eb_pad         };

    std::ofstream param_file;
    std::stringstream param_file_name;
    param_file_name << checkpointname << "/LSFactory_params";
    param_file.open(param_file_name.str());
    writeIntData(levelset_params, 4, param_file);
  }
}
