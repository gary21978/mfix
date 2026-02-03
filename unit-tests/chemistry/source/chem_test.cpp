#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#include <mfix_run_on.H>
#include <AMReX_ParticleTile.H>
#include <mfix_solvers.H>

#include <chem_test.H>

using namespace amrex;

chem_test::
~chem_test ()
{ }

chem_test::
chem_test ( )
  : m_verbose(0)
  , m_nghost(2)
  , m_newton_maxiter(500)
  , m_newton_abstol(1.e-8)
  , m_newton_reltol(1.e-8)
  , m_eb(false, maxLevel() )
  , m_level_data( maxLevel() )
{
  // Do not to iterate when creating the initial grid hierarchy
  SetIterateToFalse();

  // Use the new chopping routine which rejects cuts if
  // they don't improve the efficiency
  SetUseNewChop();

  ParmParse pp("fluid.newton_solver");

  pp.query("max_iterations", m_newton_maxiter);

  pp.query("absolute_tol", m_newton_abstol);
  pp.query("relative_tol", m_newton_reltol);

}


void
chem_test::
Init ( MFIXFluidPhase const& a_fluid,
       MFIXReactions  const& a_reactions)
{
  // Setup box array and dmap
  BoxArray ba(geom[0].Domain());
  ba.maxSize(max_grid_size[0]);
  if (ba == grids[0]) { ba = grids[0]; }

  DistributionMapping dm(ba, ParallelDescriptor::NProcs());

  finest_level = 0;
  int const lev(0);

  SetBoxArray(lev, ba);
  SetDistributionMap(lev, dm);

  // initialize default eb geometry
  m_eb.make_eb_regular(Geom());
  m_eb.make_factory(lev, Geom(lev), boxArray(lev), DistributionMap(lev));

  m_interp_data.resize(nlev());

  DualGridAuxIndexes fluid_idxs( a_fluid.solve_enthalpy(),
      a_fluid.nspecies(), a_fluid.isMixture());

  m_level_data.m_level_data[lev] =
      std::make_unique<LevelData>(a_fluid, a_reactions);

  leveldata().define(lev, m_nghost, grids[lev], dmap[lev], *(eb().factory()[lev]) );

  leveldata(lev)->initVals(0.);

  m_interp_data[lev].reset(new MultiFab(grids[lev], dmap[lev],
    fluid_idxs.count, m_nghost, MFInfo(), *(eb().factory()[lev])));
}

void
chem_test::
setup ( Vector<EBFArrayBoxFactory const*> /*a_factory*/,
        MFIXFluidPhase        const& a_fluid,
        MFIXReactions         const& /*a_reactions*/,
        MFIXInitialConditions const& a_initial_conditions)
{
  //AMREX_ALWAYS_ASSERT( m_initial_conditions.ic().size() == 1);

  if ((a_fluid.solve_enthalpy()) && a_fluid.constraint.isIdealGas() ) {
    thermo_p = a_fluid.thermodynamic_pressure();
  }

  for (int lev(0); lev<nlev(); ++lev) {

    //Print() << "Setting initial conditions for level " << lev << "\n";

    leveldata().epf_old(lev)->setVal(1., 0, 1, 0);
    leveldata().epf(lev)->setVal(1., 0, 1, 0);

    leveldata().vel_old(lev)->setVal(0., 0, 1, 0);

    Real const rho = a_initial_conditions.ic(0).fluid.get_density();
    leveldata().rho_old(lev)->setVal(rho, 0, 1, 0);

    if ( a_fluid.solve_enthalpy() ) {

      Real const Tg = a_initial_conditions.ic(0).fluid.get_temperature();
      leveldata().T_old(lev)->setVal(Tg, 0, 1, 0);

      Real hg(0.);
      {
        const auto fluid_props = a_fluid.props.data<RunOn::Host>();
        int const is_covered(0);
        if ( a_fluid.isMixture() ) {
          for (int n(0); n<a_fluid.nspecies(); ++n) {
            Real const Xn = a_initial_conditions.ic(0).fluid.get_species(n);
            hg += Xn*fluid_props.enthalpy(n, Tg, is_covered);
          }
        } else {
          hg += fluid_props.enthalpy(Tg, nullptr, is_covered);
        }
      }
      leveldata().h_old(lev)->setVal(hg, 0, 1, 0);
    }

    if ( a_fluid.isMixture() ) {

      for (int n(0); n<a_fluid.nspecies(); ++n) {
        Real const Xn = a_initial_conditions.ic(0).fluid.get_species(n);
        leveldata().X_old(lev)->setVal(Xn, n, 1, 0);
      }
    }

    // Copy old into new
    m_level_data.txfr(lev)->setVal(0.);
    m_level_data.chem_txfr(lev)->setVal(0.);


    // Copy old into new
    leveldata(lev)->resetNewWithOld();

    leveldata().txfr(lev)->setVal(0.);
    leveldata().chem_txfr(lev)->setVal(0.);

  }

}


void
chem_test::
update_interp_data ( MFIXFluidPhase const& a_fluid )
{

  DualGridAuxIndexes fluid_idxs( a_fluid.solve_enthalpy(),
      a_fluid.nspecies(), a_fluid.isMixture());

  for (int lev(0); lev<nlev(); lev++) {

    MultiFab::Copy(*m_interp_data[lev], *(leveldata().vel(lev)),
        0, fluid_idxs.vel_g, 3, m_nghost);

    MultiFab::Copy(*m_interp_data[lev], *(leveldata().epf(lev)),
        0, fluid_idxs.ep_g,  1, m_nghost);

    MultiFab::Copy(*m_interp_data[lev], *(leveldata().rho(lev)),
        0, fluid_idxs.ro_g,  1, m_nghost);

    if ( a_fluid.solve_enthalpy() ) {

      MultiFab::Copy(*m_interp_data[lev], *(leveldata().T(lev)),
          0, fluid_idxs.T_g, 1, m_nghost);

      MultiFab::Copy(*m_interp_data[lev], *(leveldata().h(lev)),
          0, fluid_idxs.h_g, 1, m_nghost);
    }

    if ( a_fluid.isMixture() ) {

      MultiFab::Copy(*m_interp_data[lev], *(leveldata().X(lev)),
          0, fluid_idxs.X_gk, a_fluid.nspecies(), m_nghost);
    }

    EB_set_covered(*m_interp_data[lev], 0, fluid_idxs.count, 0,
      std::numeric_limits<amrex::Real>::max());

    leveldata(lev)->swapOldAndNew();
  }

}
