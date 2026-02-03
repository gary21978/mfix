#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FillPatchUtil.H>

#include <mfix_filcc.H>
#include <mfix_pc.H>
#include <mfix_reporter.H>

#include <mfix_deposition_op.H>
#include <mfix_deposition_K.H>

using namespace amrex;

MFIXDepositionOp::
~MFIXDepositionOp ( )
{
  unsigned int const nlev( m_ptr.size() );
  AMREX_ASSERT( m_eps.size() == nlev );

  for (unsigned int lev(0); lev<nlev; lev++) {

    if (m_ptr[lev] != nullptr) {
      if (m_ptr[lev] != m_dst[lev]) { delete m_ptr[lev]; }
    }
    if (m_eps[lev] != nullptr) { delete m_eps[lev]; }
  }

  m_ptr.clear();
  m_eps.clear();

  m_dst.clear();
}

MFIXDepositionOp::
MFIXDepositionOp ( BCList const& a_bc_list,
                   MFIXParticleContainer* a_pc,
                   const EBFArrayBoxFactory* a_lev0_factory,
                   amrex::Vector <amrex::MultiFab*> a_dst )
  : m_bc_list(a_bc_list)
  , m_pc(a_pc)
  , m_lev0_factory(a_lev0_factory)
{

  getInputParams();

  m_dst.resize(a_dst.size());

  m_ptr.resize(m_pc->numLevels());
  m_eps.resize(m_pc->numLevels());

  // This could be removed by reworking some of the loops
  AMREX_ALWAYS_ASSERT(m_ptr.size() <= m_dst.size());

  m_on_same_grids.resize(m_pc->numLevels());

  for (int lev(0); lev<m_pc->numLevels(); ++lev) {

    // Store the pointer to the destination MultiFab.
    m_dst[lev] = a_dst[lev];

    const DistributionMapping& fluid_dm = m_dst[lev]->DistributionMap();
    const BoxArray&            fluid_ba = m_dst[lev]->boxArray();

    const DistributionMapping& particle_dm = m_pc->ParticleDistributionMap(lev);
    const BoxArray&            particle_ba = m_pc->ParticleBoxArray(lev);

    m_on_same_grids[lev] = ((fluid_dm == particle_dm) && (fluid_ba == particle_ba));

    int const ngrow(m_dst[lev]->nGrow());
    int const ncomp(m_dst[lev]->nComp());

    // If we are already working with the internal mf defined on the
    // particle_box_array, then we just work with this.
    if (lev == 0 && OnSameGrids(lev) ) {

      m_eps[lev] = new MultiFab(particle_ba, particle_dm, 1, ngrow);
      m_ptr[lev] = m_dst[lev];

    // If a_txfr is not defined on the particle_box_array, then we need
    // to make a temporary here and copy into it at the end.
    } else if (lev == 0 && !OnSameGrids(lev) ) {

      m_eps[lev] = new MultiFab(particle_ba, particle_dm, 1, ngrow);
      m_ptr[lev] = new MultiFab(particle_ba, particle_dm, ncomp, ngrow);

    } else {

      // If lev > 0 we make a temporary box array at the coarse resolution
      BoxArray crse_ba(coarsen(particle_ba, m_pc->GetParGDB()->refRatio(0)));

      m_eps[lev] = new MultiFab(crse_ba, particle_dm, 1, 1);
      m_ptr[lev] = new MultiFab(crse_ba, particle_dm, ncomp, 1);

    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( ngrow > 0,
      "MFIXDepositionOp:: Must have at least one ghost cell");

  }
}


void
MFIXDepositionOp::
CopyToDest ()
{
  Real const a_time(0.);

  // IntVect ref_ratio(this->m_gdb->refRatio(0));

  // Now interpolate from the coarse grid to define the fine grid ep-g
  Interpolater* mapper = &cell_cons_interp;
  int lo_bc[3] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
  int hi_bc[3] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
  Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

  BndryFuncArray bfunc(mfix_aux::filcc);

  for (int lev(0); lev < m_dst.size(); ++lev) {

    if (lev == 0 || !m_reduce_to_lev0) {

      // If m_dst is not defined on the particle_box_array, so we
      // copy from m_ptr into m_dst.
      if (m_dst[lev] != m_ptr[lev]) {
        m_dst[lev]->ParallelCopy(*m_ptr[lev], 0, 0, m_ptr[0]->nComp());
      }

    } else {

      PhysBCFunct<BndryFuncArray> cphysbc(m_pc->Geom(lev-1), bcs, bfunc);
      PhysBCFunct<BndryFuncArray> fphysbc(m_pc->Geom(lev  ), bcs, bfunc);

      m_dst[lev]->setVal(0.);
      InterpFromCoarseLevel(*m_dst[lev], a_time, *m_dst[lev-1], 0, 0, 1,
        m_pc->Geom(lev-1), m_pc->Geom(lev), cphysbc, 0, fphysbc, 0,
        m_pc->GetParGDB()->refRatio(0), mapper, bcs, 0);
    }
  }
}


int
MFIXDepositionOp::
getInputParams ()
{
  ParmParse pp_mfix("mfix");

  std::string deposition_scheme = "trilinear";
  pp_mfix.query("deposition_scheme", deposition_scheme);

  if (deposition_scheme.compare("trilinear") == 0) {
    m_deposition_scheme = DepositionScheme::trilinear;
  }
  else if (deposition_scheme.compare("trilinear-dpvm-square") == 0) {
    m_deposition_scheme = DepositionScheme::square_dpvm;
  }
  else if (deposition_scheme.compare("true-dpvm") == 0) {
    m_deposition_scheme = DepositionScheme::true_dpvm;
  }
  else if (deposition_scheme.compare("centroid") == 0) {
    m_deposition_scheme = DepositionScheme::centroid;
  }
  else {
    amrex::Abort("Don't know this deposition_scheme!");
  }

  m_deposition_scale_factor = 1.;
  pp_mfix.query("deposition_scale_factor", m_deposition_scale_factor);

  std::string redist_type = "MaxPack";
  pp_mfix.queryAdd("deposition_redist_type", redist_type);

  if (toLower(redist_type).compare("maxpack") == 0) {

    m_redist_type = RedistType::m_MaxPack;

    // Solids volume fractions greater than this value in cut cells
    // triggers the MaxPack redistribution algorithm.
    m_redist_max_pack = 0.6;
    pp_mfix.query("max_solids_volume_fraction", m_redist_max_pack);


  } else if (toLower(redist_type).compare("stateredist") == 0) {

    m_redist_type = RedistType::m_StateRedist;

    // Solids volume in cut cells with a geometric volume fraction below
    // m_redist_vfrac is redistributed regardless of solids volume.
    m_redist_vfrac = 0.1;
    pp_mfix.query("deposition_redist_vfrac", m_redist_vfrac);

  } else {

    reporter::Log(reporter::Error,__FILE__,__LINE__)
      << "Unknown deposition_redist_type = " << redist_type << "\n"
      << "Valid settings: MaxPack, StateRedist\n"
      << "Please correct the input deck.";
    return 1;
  }
  return 0;
}
