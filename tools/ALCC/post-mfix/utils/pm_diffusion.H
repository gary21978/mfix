#ifndef POST_MFIX_DIFFUSION_H_
#define POST_MFIX_DIFFUSION_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

// Solver for the implicit part of the diffusion equation:
//
// ( 1 - dt / rho * div ( eta grad ) ) u* = rhs
//
// Input "vel" holds the right hand side of this equation:
//
//  rhs = u + dt * ( conv + div ( eta (grad u)^T ) / rho + g - grad(p + p0) / rho )
//
// This is replaced by the solution of the matrix solve.
// Output "vel" therefore holds u*
//

class pm_diffusion {

  public:

    pm_diffusion ( int const a_max_level,
                   amrex::Vector<amrex::Geometry>            const & a_geom,
                   amrex::Vector<amrex::DistributionMapping> const & a_dmap,
                   amrex::Vector<amrex::BoxArray>            const & a_grids);

    ~pm_diffusion ();

    void reset ( amrex::Vector<amrex::MultiFab const*> const a_alpha,
                 int const a_alpha_comp, amrex::Real a_mean_diameter);

    int no_filter () { return (m_filter_type == FilterType::none) ? 1 : 0; }

    void smooth ( const amrex::Vector< amrex::MultiFab* >& a_MF,
                  int const a_srccomp, int const a_numcomp,
                  amrex::Vector<std::string> a_variables={});

  private:

    struct FilterType { enum {none, constant, variable}; };

    const int m_ngrow = 2;

    int m_max_level;

    amrex::Vector<amrex::Geometry>            const m_geom;
    amrex::Vector<amrex::DistributionMapping> const m_dmap;
    amrex::Vector<amrex::BoxArray>            const m_grids;

    int m_filter_type;

    // MLMG solver settings
    int m_verbose;
    int m_bottom_verbose;
    int m_max_iter;
    int m_bottom_max_iter;
    int m_max_fmg_iter;
    int m_linop_maxorder;

    bool m_agglomeration;
    int m_agg_grid_size;

    bool m_consolidation;
    bool m_semicoarsening;
    int m_max_coarsening_level;
    int m_max_semicoarsening_level;

    int m_pre_smooth_iter;
    int m_post_smooth_iter;
    int m_final_smooth_iter;
    int m_bottom_smooth_iter;

    amrex::Real m_rtol;
    amrex::Real m_atol;
    std::string m_bottom_solver;

    std::string m_hypre_namespace;

    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> m_bcoef;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_phi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_rhs;

    std::unique_ptr<amrex::MLABecLaplacian> m_mlabec;

  private:

    int get_nlev () { return m_max_level + 1; }

    void assert_has_mean_diameter (std::string const a_file,
                                   int const a_line,
                                   amrex::Real a_mean_diameter);

};
#endif
