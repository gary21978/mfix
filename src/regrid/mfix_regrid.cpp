#include <mfix.H>
#include <mfix_fluid.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

using namespace amrex;

void mfix::
Regrid ( int const a_step ) {
  if (m_regrid_int > 0 && (a_step%m_regrid_int == 0) ) {
    Print() << "Regridding at step " << a_step << '\n';
    Regrid();
  }
}

void mfix::
Regrid ()
{

  BL_PROFILE_REGION_START("mfix::Regrid()");

  if (max_level > 0) {
    regrid(/*base level=*/0, timer().time());
  }

  // The algorithms to 'greedy' grid particles was removed in MR1446

  int const loadbalance_fluid( fluid.solve() && m_loadbalance->UpdateFluid() );
  int const loadbalance_solids( (m_dem.solve() || m_pic.solve()) );

  { int const lev(0);

    DistributionMapping new_dm = m_loadbalance->MakeDistMapping(lev);

    if ( loadbalance_solids ) {

      if ( m_verbose > 0 ) {
        Real load_eff = pc->particleImbalance();
        Print() << "particle load efficiency before regridding " << load_eff << '\n';
      }

      pc->Regrid(new_dm, pc->ParticleBoxArray(lev), lev);

      if (sort_particle_int > 0) {
        pc->SortParticlesByBin(particle_sorting_bin);
      }

      m_eb->regrid_levelset_array(lev, Geom(), pc);

      if ( m_verbose > 0 ) {
        Real load_eff = pc->particleImbalance();
        Print() << "particle load efficiency after regridding " << load_eff << '\n';
      }
    }


    if ( loadbalance_fluid ) {

      if ( m_verbose > 0 ) { Print() << "Load balancing fluid\n"; }

      // If DualGrid and the fluid is updated, use the cell count.
      if ( m_loadbalance->DualGrid() ) {

        if ( m_verbose > 0 ) { Print() << "Making new dmap from grids\n"; }

        new_dm = m_loadbalance->
            MakeDistMappingFromGrids(grids[lev], dmap[lev]);
      }

      if ( new_dm != dmap[lev] ) {

        SetDistributionMap(lev, new_dm);

        RemakeLevel(lev, timer().time(), grids[lev], dmap[lev]);
      }
    }

    if (m_loadbalance->DualGrid() || !fluid.solve()) {

      m_loadbalance->ResetWeights(lev, pc->ParticleBoxArray(lev),
          pc->ParticleDistributionMap(lev));

    } else { AMREX_ASSERT( m_loadbalance->SingleGrid() );

      m_loadbalance->ResetWeights(lev, grids[lev], dmap[lev]);

    }

  }


  // Recompute fab areas after regrid. We only recompute the areas
  // if solving particles because the fluid doesn't use the individual
  // "fab bc areas," and the total bc area will not change.
  if ( has_particles() ) {

    for (int lev(0); lev <= pc->finestLevel(); lev++) {
      bcs().calc_bc_areas(lev, pc->ParticleBoxArray(lev),
          pc->ParticleDistributionMap(lev),
          eb()->particle_factory()[lev].get());
    }
  }

  BL_PROFILE_REGION_STOP("mfix::Regrid()");
}


//! Remake an existing level using provided BoxArray and
//! DistributionMapping and fill with existing fine and coarse data.
//! Called by AmrCore::regrid.
void mfix::
RemakeLevel ( int a_lev, Real a_time, BoxArray const& a_grids,
              DistributionMapping const& a_dmap)
{
  BL_PROFILE("mfix::RemakeLevel()");

  if (m_verbose > 0) { Print() << "Remaking level " << a_lev << '\n'; }

  eb()->update_factory(a_lev, geom[a_lev], a_grids, a_dmap);

  std::unique_ptr<LevelData> new_level( new LevelData(fluid, reactions));
  new_level->define(nghost_state(), a_grids, a_dmap, *(m_eb->factory()[a_lev]) );
  new_level->initVals();

  int const nghost(0);

  // Volume fraction
  bcs().fillpatch(a_lev, a_time, BCFillVar::epf,
      leveldata().epf(), new_level->m_epf.get(), nghost );

  // Density
  bcs().fillpatch(a_lev, a_time, BCFillVar::rho,
      leveldata().rho(), new_level->m_rho.get(), nghost );

  // Species
  if (leveldata(a_lev)->has_species()) {
    bcs().fillpatch(a_lev, a_time,  BCFillVar::X,
        leveldata().X(), new_level->m_X.get(), nghost );
  }

  // Enthalpy
  if (leveldata(a_lev)->has_enthalpy()) {
    bcs().fillpatch(a_lev, a_time, BCFillVar::h,
        leveldata().h(), new_level->m_h.get(), nghost );
  }

  // Temperature
  if (leveldata(a_lev)->has_enthalpy()) {
    bcs().fillpatch(a_lev, a_time,  BCFillVar::T,
        leveldata().T(), new_level->m_T.get(), nghost );
  }

  // Tracers
  if (leveldata(a_lev)->has_tracer()) {
    bcs().fillpatch(a_lev, a_time,  BCFillVar::tracer,
        leveldata().tracer(), new_level->m_tracer.get(), nghost );
  }

  // Velocity and pressure gradient
  bcs().fillpatch(a_lev, a_time, BCFillVar::vel,
      leveldata().vel(), new_level->m_vel.get(), nghost );

  // Velocity and pressure gradient
  bcs().fillpatch(a_lev, a_time, BCFillVar::none,
      leveldata().grad_p(), new_level->m_grad_p.get(), nghost );

  // Replace the old level with the new one we just created
  leveldata().move( a_lev, std::move(new_level) );

  // Copy data form current into old.
  leveldata(a_lev)->resetOldWithNew();

  // This call resets both the MAC, nodal and the diffusion solvers
  clear_solvers();
}


void mfix::
ClearLevel (int a_lev)
{
  leveldata().clear( a_lev );
  eb()->clear_factory( a_lev );

  // This call resets both the MAC, nodal and the diffusion solvers
  clear_solvers();
}


//! Remake an existing level using provided BoxArray and
//! DistributionMapping and fill with existing fine and coarse data.
//! Called by AmrCore::regrid.
void mfix::
MakeNewLevelFromCoarse ( int a_lev, Real a_time, BoxArray const& a_grids,
                         DistributionMapping const& a_dmap)
{
  BL_PROFILE("mfix::MakeNewLevelFromCoarse()");

  if (m_verbose > 0) { Print() << "Make new level " << a_lev << " from coarse\n"; }

  eb()->make_factory(a_lev, geom[a_lev], a_grids, a_dmap);

  std::unique_ptr<LevelData> new_level( new LevelData(fluid, reactions));
  new_level->define(nghost_state(), a_grids, a_dmap, *(m_eb->factory()[a_lev]) );
  new_level->initVals();

  int const nghost(0);

  // Volume fraction
  bcs().fillcoarsepatch(a_lev, a_time, BCFillVar::epf,
      leveldata().epf(), new_level->m_epf.get(), nghost);

  // Density
  bcs().fillcoarsepatch(a_lev, a_time, BCFillVar::rho,
      leveldata().rho(), new_level->m_rho.get(), nghost);

  // Species
  if (leveldata(a_lev)->has_species()) {
    bcs().fillcoarsepatch(a_lev, a_time,  BCFillVar::X,
        leveldata().X(), new_level->m_X.get(), nghost);
  }

  // Enthalpy
  if (leveldata(a_lev)->has_enthalpy()) {
    bcs().fillcoarsepatch(a_lev, a_time, BCFillVar::h,
        leveldata().h(), new_level->m_h.get(), nghost );
  }

  // Temperature
  if (leveldata(a_lev)->has_enthalpy()) {
    bcs().fillcoarsepatch(a_lev, a_time,  BCFillVar::T,
        leveldata().T(), new_level->m_T.get(), nghost );
  }

  // Tracers
  if (leveldata(a_lev)->has_tracer()) {
    bcs().fillcoarsepatch(a_lev, a_time,  BCFillVar::tracer,
        leveldata().tracer(), new_level->m_tracer.get(), nghost );
  }

  // Velocity and pressure gradient
  bcs().fillcoarsepatch(a_lev, a_time, BCFillVar::vel,
      leveldata().vel(), new_level->m_vel.get(), nghost );

  // Velocity and pressure gradient
  bcs().fillcoarsepatch(a_lev, a_time, BCFillVar::none,
      leveldata().grad_p(), new_level->m_grad_p.get(), nghost );

  // This call resets both the MAC, nodal and the diffusion solvers
  clear_solvers();
}


void mfix::
clear_solvers ()
{
  // Reset macproj solver
  if ( macproj ) { macproj.reset(); }

  // Reset diffusion solvers
  if ( m_diffOpVel      ) { m_diffOpVel.reset();      }
  if ( m_diffOpSpecies  ) { m_diffOpSpecies.reset();  }
  if ( m_diffOpEnergy   ) { m_diffOpEnergy.reset();   }
  if ( m_diffOpTracer   ) { m_diffOpTracer.reset();   }
  if ( m_diffOpVoidFrac ) { m_diffOpVoidFrac.reset(); }
  if ( m_diffOpTxfr     ) { m_diffOpTxfr.reset();     }
}
