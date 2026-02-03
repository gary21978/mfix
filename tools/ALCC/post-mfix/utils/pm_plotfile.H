#ifndef POST_MFIX_PLOTFILE_
#define POST_MFIX_PLOTFILE_

#include <AMReX_AmrMesh.H>
#include <AMReX_Geometry.H>
#include <AMReX_PlotFileUtil.H>

class pm_plotfile : protected amrex::PlotFileData
{
  public:

    pm_plotfile ( std::string a_pf_name );

    std::string get_name () noexcept { return m_pf_name; }

    amrex::AmrInfo get_amr_info () noexcept { return m_amr_info; }
    amrex::Geometry get_level_0_geom () noexcept { return m_level_0_geom; }

    int get_nlev () const noexcept { return finestLevel() + 1; }

    int get_step () const noexcept { return levelStep(0); }

    amrex::Real get_time () const noexcept { return time(); }

    const amrex::BoxArray& get_ba (int a_lev) const noexcept
    { return boxArray(a_lev); }

    const amrex::DistributionMapping& get_dm (int a_lev) const noexcept
    { return DistributionMap(a_lev); }

    int refRatio (int a_lev) const noexcept { return refRatio(a_lev); }

    // Methods of accessing basic plot file grid data.
    //-------------------------------------------------------------------------

    // Number of fluid (grid) components
    int fluid_comps () const noexcept { return nComp(); }

    // Names of fluid variables
    const amrex::Vector<std::string>& fluid_variables () const noexcept
    { return varNames(); }

    // Access to level data in the plot file
    amrex::MultiFab get_mf_data (int a_lev) noexcept { return get(a_lev); }


    // Methods of accessing basic plot file particle data.
    //-------------------------------------------------------------------------

    int particle_count () noexcept { return m_nparticles; }

    int particle_real_comps () noexcept { return m_pf_size_r_extra; }

    const amrex::Vector<std::string>& particle_real_variables () const noexcept
    { return m_pf_variables_r; }

    int particle_int_comps () noexcept { return m_pf_size_i; }

    const amrex::Vector<std::string>& particle_int_variables () const noexcept
    { return m_pf_variables_i; }


  private:

    std::string m_pf_name;

    amrex::Geometry m_level_0_geom;

    amrex::AmrInfo m_amr_info;

    // dimension
    int m_pf_dim;

    int m_pf_size_r_extra; // number of reals, excluding position

    // Variable names (reals)
    amrex::Vector<std::string> m_pf_variables_r;

    int m_pf_size_i;       // number of ints, excluding id and cpu

    // Variable names (ints)
    amrex::Vector<std::string> m_pf_variables_i;

    bool m_is_checkpoint;

    int m_nparticles;

    // Operator to read in some header information.
    friend std::istream& operator>> (std::istream& a_stream, pm_plotfile& header);
};
#endif
