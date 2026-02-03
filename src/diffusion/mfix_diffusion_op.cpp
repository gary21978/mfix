#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_utils.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_EB_Redistribution.H>

#include <mfix_diffusion_op.H>
#include <mfix_utils_harmonic_avg_K.H>

#include <mfix_run_on.H>

using namespace amrex;

MFIXDiffusionOp::
MFIXDiffusionOp ( int const a_nlev,
                  Vector<Geometry> const& a_geom,
                  int const a_ncomp,
                  MFIXFluidPhase& a_fluid,
                  MFIXBoundaryConditions& a_bcs,
                  Vector<BCRec> const& a_bcrec)
  : m_nlev(a_nlev)
  , m_geom(a_geom)
  , m_ncomp(a_ncomp)
  , m_fluid(a_fluid)
  , m_bcs(a_bcs)
  , m_bcrec(a_bcrec)
{
}


void MFIXDiffusionOp::
readParameters ()
{
  ParmParse pp("diffusion");

  pp.query("verbose_solver",          m_verbose);

  pp.query("verbose",                 m_mg_verbose);
  pp.query("bottom_verbose",          m_mg_bottom_verbose);

  pp.query("maxiter",                 m_mg_maxiter);
  pp.query("bottom_maxiter",          m_mg_bottom_maxiter);
  pp.query("mg_max_coarsening_level", m_mg_max_coarsening_level);
  pp.query("rtol",                    m_mg_rtol);
  pp.query("atol",                    m_mg_atol);

  pp.query("mg_max_fmg_iter",         m_mg_max_fmg_iter);
  pp.query("agg_grid_size",           m_mg_agg_grid_size);

  pp.query("bottom_solver",           m_bottom_solver_type);
  pp.query("hypre_namespace",         m_hypre_namespace);
  pp.query("hypre_interface",         m_hypre_interface);
}


void MFIXDiffusionOp::
define ( Vector<BoxArray>                   const& a_grids,
         Vector<DistributionMapping>        const& a_dmap,
         Vector< const EBFArrayBoxFactory* >const& a_ebfactory,
         int const a_matrix_ncomp )
{

  m_b.resize(m_nlev);
  m_phi.resize(m_nlev);
  m_rhs.resize(m_nlev);

  for(int lev(0); lev<m_nlev; lev++) {

    for(int idim(0); idim < AMREX_SPACEDIM; idim++) {

      BoxArray edge_ba = a_grids[lev];
      edge_ba.surroundingNodes(idim);

      m_b[lev][idim].reset( new MultiFab(edge_ba, a_dmap[lev],
          a_matrix_ncomp, 0, MFInfo(), *a_ebfactory[lev]));

      m_b[lev][idim]->setVal(0.);
    }

    m_phi[lev].reset( new MultiFab(a_grids[lev], a_dmap[lev],
        m_ncomp, 1, MFInfo(), *a_ebfactory[lev]));

    m_phi[lev]->setVal(0.);

    // No ghost cells needed for rhs
    m_rhs[lev].reset( new MultiFab(a_grids[lev], a_dmap[lev],
        m_ncomp, 0, MFInfo(), *a_ebfactory[lev]));

    m_rhs[lev]->setVal(0.);

  }
}


//
// Set the user-supplied settings for the MLMG solver
// (this must be done every time step, since MLMG is created after updating matrix
//
void MFIXDiffusionOp::
setSolverSettings (MLMG& a_solver)
{
    // The default bottom solver is BiCG
  if(m_bottom_solver_type == "smoother") {

    a_solver.setBottomSolver(MLMG::BottomSolver::smoother);

  } else if(m_bottom_solver_type == "hypre") {
#ifdef AMREX_USE_HYPRE
    a_solver.setBottomSolver(MLMG::BottomSolver::hypre);
    a_solver.setHypreOptionsNamespace(m_hypre_namespace);
    if (m_hypre_interface == "ij") {
       a_solver.setHypreInterface(amrex::Hypre::Interface::ij);
    } else if (m_hypre_interface == "semi_structured") {
       a_solver.setHypreInterface(amrex::Hypre::Interface::semi_structed);
    } else if (m_hypre_interface == "structured") {
       a_solver.setHypreInterface(amrex::Hypre::Interface::structed);
    } else {
     amrex::Abort(
           "Invalid hypre interface. Valid options: ij semi_structured "
           "structured");
    }
#else
    amrex::Abort("mfix was not built with hypre support");
#endif
  }

  // Maximum iterations for MultiGrid / ConjugateGradients
  a_solver.setMaxIter(m_mg_maxiter);
  a_solver.setMaxFmgIter(m_mg_max_fmg_iter);
  a_solver.setBottomMaxIter(m_mg_bottom_maxiter);

  // Verbosity for MultiGrid / ConjugateGradients
  a_solver.setVerbose(m_mg_verbose);
  a_solver.setBottomVerbose(m_mg_bottom_verbose);

  // This ensures that ghost cells of phi are correctly filled when
  // returned from the solver
  a_solver.setFinalFillBC(true);
}


void MFIXDiffusionOp::
setDiffCoeff( Vector< MultiFab const*> const& a_b_cc )
{
  for (auto& lev_mf : m_b_cc) {
    if (lev_mf != nullptr ) { lev_mf.reset(); }
  }
  AMREX_ASSERT( m_nlev == a_b_cc.size() );

  // We could loosen this restriction, but will need to make a similar
  // change to the matrix definition.
  AMREX_ASSERT( a_b_cc[0]->nComp() == 1 );

  m_b_cc.resize(m_nlev);
  for (int lev(0); lev<m_nlev; ++lev) {

    m_b_cc[lev].reset( new MultiFab(m_rhs[lev]->boxArray(),
        m_rhs[lev]->DistributionMap(), 1, 1, MFInfo(),
        m_rhs[lev]->Factory()) );

    MultiFab::Copy(*m_b_cc[lev], *a_b_cc[lev], 0, 0, 1, m_b_cc[lev]->nGrow());
  }
  m_const_coeff = 0;
  m_b_defined   = 1;
}


void MFIXDiffusionOp::
setDiffCoeff( Real a_b_const ) {
  Vector<Real> b_const = { a_b_const };
  setDiffCoeff(b_const);
}


void MFIXDiffusionOp::
setDiffCoeff( Vector<Real> a_b_const ) {

  // This could be relaxed with a few additional changes.
  AMREX_ASSERT( a_b_const.size() == 1 );
  //AMREX_ASSERT( a_b_const.size() == m_ncomp);

  m_b_const.resize(m_ncomp);
  int const nval( a_b_const.size() );
  for (int n(0); n<m_ncomp; ++n) {
    if (nval == 1) { m_b_const[n] = a_b_const[0]; }
    else { m_b_const[n] = a_b_const[n]; }
  }
  m_b_defined   = 1;
}


void MFIXDiffusionOp::
setEBDirichlet ()
{
  m_set_eb_dirichlet = 1;

  m_b_eb.resize(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    m_b_eb[lev].reset( new MultiFab( m_rhs[lev]->boxArray(),
        m_rhs[lev]->DistributionMap(), 1, 1, MFInfo(),
        m_rhs[lev]->Factory()) );

    m_b_eb[lev]->setVal(0.);
  }
}


//
// Set the BCs for face-centroid-based velocity components only
//
void
MFIXDiffusionOp::
avgDiffCoeffToFaces ( int const a_lev )
{
  BL_PROFILE("MFIXDiffusionOp::avgDiffCoeffToFaces()");

  const auto& factory =
      dynamic_cast<EBFArrayBoxFactory const&>(m_b_cc[a_lev]->Factory());

  const auto& flags = factory.getMultiEBCellFlagFab();
  const auto& vfrac = factory.getVolFrac();
  const auto& area  = factory.getAreaFrac();
  const auto& ccent = factory.getCentroid();

  Box domain(m_geom[a_lev].Domain());

  for (MFIter mfi(*m_b_cc[a_lev], false); mfi.isValid(); ++mfi) {

    Box const& bx = mfi.tilebox();

    const Box& xbx = mfi.nodaltilebox(0);
    const Box& ybx = mfi.nodaltilebox(1);
    const Box& zbx = mfi.nodaltilebox(2);

    const auto fabtyp = flags[mfi].getType(amrex::grow(bx,0));
    const auto fabtyp_ghost = flags[mfi].getType(amrex::grow(bx,1));

    int const ncomp(m_ncomp);

    int const xp(m_geom[a_lev].isPeriodic(0));
    int const yp(m_geom[a_lev].isPeriodic(1));
    int const zp(m_geom[a_lev].isPeriodic(2));

    if (fabtyp != FabType::covered) {

      Array4<Real const> const& b_cc = m_b_cc[a_lev]->const_array(mfi);

      Array4<Real> const& b_x = m_b[a_lev][0]->array(mfi);
      Array4<Real> const& b_y = m_b[a_lev][1]->array(mfi);
      Array4<Real> const& b_z = m_b[a_lev][2]->array(mfi);

      if (fabtyp_ghost == FabType::regular ) {

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
        (xbx, txbx, { harmonic_avg_cc2face_x(txbx, ncomp, b_x, b_cc, xp, domain); },
         ybx, tybx, { harmonic_avg_cc2face_y(tybx, ncomp, b_y, b_cc, yp, domain); },
         zbx, tzbx, { harmonic_avg_cc2face_z(tzbx, ncomp, b_z, b_cc, zp, domain); });

      } else {

        Array4<Real const> const& apx = area[0]->const_array(mfi);
        Array4<Real const> const& apy = area[1]->const_array(mfi);
        Array4<Real const> const& apz = area[2]->const_array(mfi);

        Array4<Real const> const& cvol = vfrac.const_array(mfi);
        Array4<Real const> const& cct  = ccent.const_array(mfi);

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
        (xbx, txbx, { eb_harmonic_avg_cc2face_x(txbx, ncomp, b_x, b_cc, apx, cvol, cct, xp, domain); },
         ybx, tybx, { eb_harmonic_avg_cc2face_y(tybx, ncomp, b_y, b_cc, apy, cvol, cct, yp, domain); },
         zbx, tzbx, { eb_harmonic_avg_cc2face_z(tzbx, ncomp, b_z, b_cc, apz, cvol, cct, zp, domain); });

      }
    }

  } // MFIter

  m_b[a_lev][0]->FillBoundary(m_geom[a_lev].periodicity());
  m_b[a_lev][1]->FillBoundary(m_geom[a_lev].periodicity());
  m_b[a_lev][2]->FillBoundary(m_geom[a_lev].periodicity());

}
