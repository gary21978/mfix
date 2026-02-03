#ifndef DIFFUSION_EQUATION_H_
#define DIFFUSION_EQUATION_H_

#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLEBTensorOp.H>
#include <AMReX_EBFArrayBox.H>

#include <mfix_bc.H>
#include <mfix_fluid.H>
#include <mfix_avg_pc_parms.H>
#include <mfix_reporter.H>

//
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
// Note: we actually solve the above equation multiplied by the density ro.
//
class MFIXDiffusionOp
{
  public:

    // Constructor, with initialisation
    MFIXDiffusionOp ( int const a_nlev,
                      amrex::Vector<amrex::Geometry> const& a_geom, int const a_ncomp,
                      MFIXFluidPhase& a_fluid, MFIXBoundaryConditions& a_bcs,
                      amrex::Vector<amrex::BCRec> const& a_bcrec);

    // Destructor
    virtual ~MFIXDiffusionOp () = default;

    virtual void setDiffCoeff( amrex::Real a_b_const );

    virtual void setDiffCoeff( amrex::Vector<amrex::Real> a_b_const );

    virtual void setDiffCoeff( amrex::Vector< amrex::MultiFab const*> const& a_b_cc );

    /// Take the harmonic average of m_b_cc to compute value at faces m_b
    void avgDiffCoeffToFaces (int const a_lev);

    amrex::Vector<amrex::MultiFab const*> getDiffCoeff ( ) {
      AMREX_ASSERT( m_b_defined == 1 );
      return GetVecOfConstPtrs(m_b_cc);
    }

    void setEBDirichlet ();

  protected:

    // Read input from ParmParse (solver settings)
    void readParameters ();

    void define ( amrex::Vector<amrex::BoxArray>                    const& a_grids,
                  amrex::Vector<amrex::DistributionMapping>         const& a_dmap,
                  amrex::Vector< const amrex::EBFArrayBoxFactory* > const& a_ebfactory,
                  int const a_matrix_ncomp );

    /// Set user-supplied solver settings (must be done every time step)
    void setSolverSettings (amrex::MLMG& solver);

    int m_verbose = 0;

    int const m_nlev;

    amrex::Vector<amrex::Geometry> m_geom;

    int const m_ncomp;

    MFIXFluidPhase& m_fluid;

    MFIXBoundaryConditions& m_bcs;

    amrex::Vector<amrex::BCRec> const& m_bcrec;

    // Internal data used in the matrix solve
    //
    // MultiLevel EmbeddedBoundary cell-centered Laplacian:
    //
    // ( alpha a - beta div ( b grad ) ) phi = rhs
    //

    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> m_b;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_rhs;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_phi;

    // Options to control MLMG behavior
    int m_mg_verbose = 0;
    int m_mg_bottom_verbose = 0;
    int m_mg_maxiter = 100;
    int m_mg_bottom_maxiter = 100;
    int m_mg_max_fmg_iter = 0;
    int m_mg_max_coarsening_level = 100;
    int m_mg_agg_grid_size = -1;

    amrex::Real m_mg_rtol = 1.0e-11;
    amrex::Real m_mg_atol = 1.0e-14;

    std::string m_bottom_solver_type = "bicgstab";
    std::string m_hypre_namespace = "hypre";
    std::string m_hypre_interface = "ij";

    /// Flag to set the EB Dirichlet value
    int m_set_eb_dirichlet = 0;

    int m_b_defined   = 0;
    int m_const_coeff = 1;

    amrex::Vector< amrex::Real > m_b_const;

    /// cell-centered diffusion coefficient
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_b_cc;

    /// EB Dirichlet diffusion coefficient
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_b_eb;

};


class MFIXDiffOpVelocity : public MFIXDiffusionOp
{
  public:

    // Constructor, with initialisation
    MFIXDiffOpVelocity ( int const a_nlev,
                         amrex::Vector<amrex::Geometry> const& a_geom,
                         amrex::Vector<amrex::BoxArray> const& a_grids,
                         amrex::Vector<amrex::DistributionMapping> const& a_dmap,
                         amrex::Vector< const amrex::EBFArrayBoxFactory* >const& a_ebfactory,
                         MFIXFluidPhase& a_fluid, MFIXAvgParticleParms& a_avg_pc_parms,
                         MFIXBoundaryConditions& a_bcs,
                         amrex::Vector<amrex::BCRec> const& a_bcrec);

    ~MFIXDiffOpVelocity () = default;

    void setDiffCoeff( amrex::Real /*a_b_const*/ ) override
    { reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "MFIXDiffOpVelocity class does not support the setDiffCoeff method.";
    }

    void setDiffCoeff( amrex::Vector<amrex::Real> /*a_b_const*/ ) override
    { reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "MFIXDiffOpVelocity class does not support the setDiffCoeff method.";
    }

    void setDiffCoeff( amrex::Vector< amrex::MultiFab const*> const& /*a_b_cc*/ ) override
    { reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "MFIXDiffOpVelocity class does not support the setDiffCoeff method.";
    }

    void setDiffCoeff (amrex::Vector<amrex::MultiFab const* >const& a_vel,
                       amrex::Vector<amrex::MultiFab const* >const& a_epf,
                       amrex::Vector<amrex::MultiFab const* >const& a_rho,
                       amrex::Vector<amrex::MultiFab const* >const& a_Tf,
                       amrex::Vector<amrex::MultiFab const* >const& a_avg_particle_data,
                       amrex::Vector< amrex::MultiFab const* >const& a_Xfk,
                       int const a_include_eddy_viscosity=1);

    // Explicit tensor viscous term
    void computeDivTau (amrex::Vector< amrex::MultiFab      * > const& a_divtau,
                        amrex::Vector< amrex::MultiFab      * > const& a_vel,
                        amrex::Vector< amrex::MultiFab const* > const& a_epf,
                        amrex::Vector< amrex::MultiFab const* > const& a_rho,
                        amrex::Vector< amrex::MultiFab const* > const& a_Tf,
                        amrex::Vector< amrex::MultiFab const* > const& a_avg_particle_data,
                        amrex::Vector< amrex::MultiFab const* > const& X_gk_in,
                        amrex::Vector< amrex::MultiFab const* > const& a_eb_inflow,
                        int const a_include_eddy_viscosity=1 );

    void solve ( amrex::Real const a_dt, int const a_drag_includes_divtau,
                 amrex::Vector< amrex::MultiFab      * >const& a_vel,
                 amrex::Vector< amrex::MultiFab const* >const& a_epf,
                 amrex::Vector< amrex::MultiFab const* >const& a_rho,
                 amrex::Vector< amrex::MultiFab const* >const& a_Tf,
                 amrex::Vector< amrex::MultiFab const* >const& a_avg_particle_data,
                 amrex::Vector< amrex::MultiFab const* >const& a_Xfk,
                 amrex::Real const a_drag_dt,
                 amrex::Vector< amrex::MultiFab const* >const& a_S_p,
                 amrex::Vector< amrex::MultiFab const* >const& a_S_c,
                 amrex::Vector< amrex::MultiFab const* >const& a_eb_inflow);

  private:

    MFIXAvgParticleParms& m_avg_pc_parms;

    std::unique_ptr<amrex::MLEBTensorOp> m_matrix;

};


class MFIXDiffOpTracer : public MFIXDiffusionOp
{
  public:

    // Constructor, with initialisation
    MFIXDiffOpTracer ( int const a_nlev,
                       amrex::Vector<amrex::Geometry> const& a_geom,
                       amrex::Vector<amrex::BoxArray> const& a_grids,
                       amrex::Vector<amrex::DistributionMapping> const& a_dmap,
                       amrex::Vector< const amrex::EBFArrayBoxFactory* >const& a_ebfactory,
                       int const a_ncomp, int const a_matrix_ncomp,
                       MFIXFluidPhase& a_fluid, MFIXBoundaryConditions& a_bcs,
                       amrex::Vector<amrex::BCRec> const& a_bcrec);

    ~MFIXDiffOpTracer () = default;

    // Explicit tensor viscous term
    void computeLap ( amrex::Vector< amrex::MultiFab*      > const& a_lapT,
                      amrex::Vector< amrex::MultiFab*      > const& a_scalar,
                      amrex::Vector< amrex::MultiFab const*> const& a_epf,
                      amrex::Vector< amrex::MultiFab const*> const& a_rho);

    // Implicit diffusion solve
    void solve ( amrex::Vector< amrex::MultiFab      * > const& a_scalar,
                 amrex::Vector< amrex::MultiFab const* > const& a_epf,
                 amrex::Vector< amrex::MultiFab const* > const& a_rho,
                 amrex::Real const a_dt);

    std::unique_ptr<amrex::MLEBABecLap> m_matrix;


};


class MFIXDiffOpSpecies : public MFIXDiffusionOp
{
  public:

    // Constructor, with initialisation
    MFIXDiffOpSpecies ( int const a_nlev,
                        amrex::Vector<amrex::Geometry> const& a_geom,
                        amrex::Vector<amrex::BoxArray> const& a_grids,
                        amrex::Vector<amrex::DistributionMapping> const& a_dmap,
                        amrex::Vector< const amrex::EBFArrayBoxFactory* >const& a_ebfactory,
                        int const a_ncomp, MFIXFluidPhase& a_fluid, MFIXBoundaryConditions& a_bcs,
                        amrex::Vector<amrex::BCRec> const& a_bcrec);

    ~MFIXDiffOpSpecies () = default;

    void setDiffCoeff( amrex::Real /*a_b_const*/ ) override
    { reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "MFIXDiffOpSpecies class does not support the setDiffCoeff method.";
    }

    void setDiffCoeff( amrex::Vector<amrex::Real> /*a_b_const*/ ) override
    { reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "MFIXDiffOpSpecies class does not support the setDiffCoeff method.";
    }

    void setDiffCoeff( amrex::Vector< amrex::MultiFab const*> const& /*a_b_cc*/ ) override
    { reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "MFIXDiffOpSpecies class does not support the setDiffCoeff method.";
    }

    void setDiffCoeff( int const a_lev,
                       amrex::Vector< amrex::MultiFab const*> const& a_epf,
                       amrex::Vector< amrex::MultiFab const*> const& a_rho);

    // Explicit tensor species mass fractions term
    void computeFlux (const amrex::Vector< amrex::Array< amrex::MultiFab*, AMREX_SPACEDIM> >& J_gk,
                      const amrex::Vector< amrex::MultiFab*      >& X_gk_in,
                      const amrex::Vector< amrex::MultiFab const*>& ro_g_in,
                      const amrex::Vector< amrex::MultiFab const*>& ep_g_in);

    // Computes the divergence of the species mass flux
    void computeDivJ (const amrex::Vector< amrex::MultiFab*      >& divJ_out,
                      const amrex::Vector< amrex::Array< amrex::MultiFab*, AMREX_SPACEDIM> >& J_gk);

    // Species mass fraction diffusion
    void solve ( amrex::Vector< amrex::MultiFab      * > const& a_X_gk,
                 amrex::Vector< amrex::Array< amrex::MultiFab*, AMREX_SPACEDIM>> const& J_gk,
                 amrex::Vector< amrex::MultiFab const* > const& a_epf,
                 amrex::Vector< amrex::MultiFab const* > const& a_rho,
                 amrex::Real const a_dt);

    std::unique_ptr<amrex::MLEBABecLap> m_matrix;

};


class MFIXDiffOpEnergy : public MFIXDiffusionOp
{
  public:

    // Constructor, with initialisation
    MFIXDiffOpEnergy ( int const a_nlev,
                       amrex::Vector<amrex::Geometry> const& a_geom,
                       amrex::Vector<amrex::BoxArray> const& a_grids,
                       amrex::Vector<amrex::DistributionMapping> const& a_dmap,
                       amrex::Vector< const amrex::EBFArrayBoxFactory* >const& a_ebfactory,
                       MFIXFluidPhase& a_fluid, MFIXBoundaryConditions& a_bcs,
                       amrex::Vector<amrex::BCRec> const& a_bcrec);

    ~MFIXDiffOpEnergy () = default;

    void setDiffCoeff( amrex::Real /*a_b_const*/ ) override
    { reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "MFIXDiffOpEnergy class does not support the setDiffCoeff method.";
    }

    void setDiffCoeff( amrex::Vector<amrex::Real> /*a_b_const*/ ) override
    { reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "MFIXDiffOpEnergy class does not support the setDiffCoeff method.";
    }

    void setDiffCoeff( amrex::Vector< amrex::MultiFab const*> const& /*a_b_cc*/ ) override
    { reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "MFIXDiffOpEnergy class does not support the setDiffCoeff method.";
    }

    void setDiffCoeff( int const a_lev,
                       amrex::Vector< amrex::MultiFab const*> const& a_epf,
                       amrex::Vector< amrex::MultiFab const*> const& a_Tf );

    void setEBDiffCoeff( int const a_lev,
                         amrex::Vector< amrex::MultiFab const*> const& a_Teb );

    // Explicit tensor temperature term
    void computeLapT (amrex::Vector< amrex::MultiFab*      > const& a_lapT,
                      amrex::Vector< amrex::MultiFab*      > const& a_Tf,
                      amrex::Vector< amrex::MultiFab const*> const& a_epf,
                      amrex::Vector< amrex::MultiFab const*> const& a_eb_Tf);

    // Temperature diffusion
    void solve ( amrex::Vector< amrex::MultiFab      * > const& a_Tf,
                 amrex::Vector< amrex::MultiFab      * > const& a_hf,
                 amrex::Vector< amrex::MultiFab const* > const& a_epf,
                 amrex::Vector< amrex::MultiFab const* > const& a_rho,
                 amrex::Vector< amrex::MultiFab const* > const& a_X_f,
                 amrex::Vector< amrex::MultiFab const*> const& a_eb_Tf,
                 amrex::Real const a_dt,
                 amrex::Real const a_abstol,
                 amrex::Real const a_reltol,
                 int const a_maxiter);

    // Computes the divergence of tenthalpy flux due to species mass flux
    void computeDivhJ (amrex::Vector< amrex::MultiFab*      > const& a_divhJ,
                       amrex::Vector< amrex::Array< amrex::MultiFab*, AMREX_SPACEDIM> > const& a_h_fc,
                       amrex::Vector< amrex::Array< amrex::MultiFab*, AMREX_SPACEDIM> > const& a_J_k,
                       amrex::Vector< amrex::MultiFab const*> const& a_Tf,
                       const int a_update_enthalpies = 0);

    std::unique_ptr<amrex::MLEBABecLap> m_matrix;
};


class MFIXDiffOpScalar : public MFIXDiffusionOp
{
  public:

    // Constructor, with initialisation
    MFIXDiffOpScalar ( int const a_nlev,
                       amrex::Vector<amrex::Geometry> const& a_geom,
                       amrex::Vector<amrex::BoxArray> const& a_grids,
                       amrex::Vector<amrex::DistributionMapping> const& a_dmap,
                       amrex::Vector< const amrex::EBFArrayBoxFactory* >const& a_ebfactory,
                       int const a_ncomp, int const a_matrix_ncomp,
                       MFIXFluidPhase& a_fluid, MFIXBoundaryConditions& a_bcs,
                       amrex::Vector<amrex::BCRec> const& a_bcrec);

    ~MFIXDiffOpScalar () = default;

    // Explicit tensor viscous term
    void computeLap ( amrex::Vector< amrex::MultiFab*      > const& a_lapT,
                      amrex::Vector< amrex::MultiFab*      > const& a_scalar,
                      amrex::Vector< amrex::MultiFab const*> const& a_epf,
                      amrex::Vector< amrex::MultiFab const*> const& a_rho);

    // Implicit diffusion solve
    void solve ( amrex::Vector< amrex::MultiFab      * > const& a_scalar);

    std::unique_ptr<amrex::MLEBABecLap> m_matrix;

};

#endif
