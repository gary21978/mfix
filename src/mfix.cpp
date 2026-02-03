#include <mfix.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <mfix_fluid.H>
#include <mfix_species.H>

RunType mfix::m_run_type       = RunType::Standard;

RealVect mfix::gravity {0.};
RealVect mfix::gp0     {0.};

// Destructor
mfix::~mfix ()
{
  if (m_rw != nullptr) { delete m_rw; }
  if (m_eb != nullptr) { delete m_eb; }

}

// Constructor
mfix::mfix ()
  : m_boundary_conditions(maxLevel(), Geom(), refRatio(), bc_list, m_embedded_boundaries)
  , bc_list(maxLevel() + 1)
  , m_level_data( maxLevel() )
{
  // NOTE: Geometry on all levels has just been defined in the AmrCore
  // constructor. No valid BoxArray and DistributionMapping have been defined.
  // But the arrays for them have been resized.
  Print() << "Number of levels: " << maxLevel() + 1 << '\n';

  // We have to check on both mfix and amr prefix for restart flags
  // to support the CCSE regtest.
  if ( ParmParse("mfix").contains("restart") ||
       ParmParse("amr").contains("restart") ) {
    m_run_type = RunType::Restart;
    Print() << "RunType is restart\n";

  } else if (ParmParse("pic2dem").contains("convert")) {
    m_run_type = RunType::PIC2DEM;
    Print() << "RunType is pic2dem conversion\n";

  } else {
    m_run_type = RunType::Standard;
    Print() << "RunType is Standard\n";
  }

  m_eb = new MFIXEB(ooo_debug, maxLevel());

  m_mass_balance = std::make_shared<MFIXMassBalance>(geom, m_eb->factory(),
      leveldata(), pc, fluid, solids, m_dem, m_pic, bc_list);

  m_rw = new MFIXReadWrite(maxLevel()+1, grids, geom, pc, fluid, leveldata(), therm_p,
                           m_eb->factory(), dmap, ooo_debug, m_eb->level_sets(),
                           m_eb->levelset_refinement(), m_eb->levelset_pad(),
                           m_eb->levelset_eb_refinement(), m_eb->levelset_eb_pad(), solids,
                           m_dem, m_pic, reactions, refRatio(),
                           bc_list, m_eb->particle_factory(), m_mass_balance,
                           regions, m_boundary_conditions, bcs().get_velocity_bcrec() );
}


void
mfix::avgDown (int crse_lev, const MultiFab& S_fine, MultiFab& S_crse)
{
    BL_PROFILE("mfix::avgDown()");

    amrex::EB_average_down(S_fine, S_crse, 0, S_fine.nComp(), refRatio(crse_lev));
}
