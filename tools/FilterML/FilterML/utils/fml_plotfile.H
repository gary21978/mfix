#ifndef FILTER_ML_PLOTFILE_
#define FILTER_ML_PLOTFILE_

#include <AMReX_Geometry.H>
#include <AMReX_PlotFileUtil.H>

class fml_plotfile : public amrex::PlotFileData
{
  public:

    fml_plotfile ( std::string a_pf_name );

    std::string get_name () noexcept { return m_pf_name; }

    bool has_particles () noexcept { return m_has_particles; }

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

    bool m_has_particles;

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
    friend std::istream& operator>> (std::istream& a_stream, fml_plotfile& header);
};
#endif
