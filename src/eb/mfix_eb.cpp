#include <AMReX_EB2.H>
#include <AMReX_EB_utils.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>

#include <mfix_eb.H>
#include <mfix_bc.H>
#include <mfix_reporter.H>

#include <algorithm>

using namespace amrex;

MFIXEB::MFIXEB (bool debug,
                int max_level)
  : m_support_level(EBSupport::full)
  , m_debug(debug)
  , m_max_level(max_level)
{
  ParmParse pp("mfix");

  // Parameters used be the level-set algorithm. Refer to LSFactory (or
  // mfix.H) for more details:
  //   -> refinement: how well resolved (fine) the (level-set/EB-facet)
  //                  grid needs to be (note: a fine level-set grid means
  //                  that distances and normals are computed accurately)
  //   -> pad:        how many (refined) grid points _outside_ the
  //                  problem domain the grid extends (avoids edge cases
  //                  in physical domain)
  pp.query("levelset__refinement", m_levelset_refinement);

  // Not needed here... the role of refining EB is filled with AMR level-set
  m_levelset_eb_refinement = 1;
  // Make sure that a coarsened level-set has a level-set pad of _at least_ 2;
  m_levelset_pad = 2*m_levelset_refinement;
  // Ensure that velocity_reconstruction has enough level-set to work off:
  // (2 => EB lives on the same grid resolution as fluid)
  m_levelset_eb_pad = amrex::max(2, m_levelset_pad);

  reporter::Log(reporter::Status)
      << "Auto-generating level-set parameters:"
      << "\n  eb_refinement = " << m_levelset_eb_refinement
      << "\n  levelset_pad  = " << m_levelset_pad
      << "\n  eb_pad        = " << m_levelset_eb_pad;

  pp.query("write_eb_surface", m_write_surface);

  m_ebfactory.resize(m_max_level+1);

  m_particle_ebfactory.resize(m_max_level+1);

  m_eb_levels.resize(amrex::max(2, m_max_level+1));

  m_level_sets.resize(amrex::max(2, m_max_level+1));

}

void MFIXEB::
make_geometry ( Vector<Geometry> const& a_geom,
                std::string const& restart_file)
{
  /****************************************************************************
   *                                                                          *
   * mfix.geometry=<string> specifies the EB geometry.                        *
   * <string> can be one of: box, cylinder, hopper, generic (or blank)        *
   *                                                                          *
   ***************************************************************************/

  ParmParse pp("mfix");

  // EB geometry defined by user
  std::string geom_type = "box";
  pp.queryAdd("geometry", geom_type);

  geom_type = toLower(geom_type); // case insensitive

  // Geometry checkpoint read/write flags
  bool geom_chk_read = false;
  pp.queryAdd("geom_chk_read", geom_chk_read);

  bool geom_chk_write = false;
  pp.queryAdd("geom_chk_write", geom_chk_write);

  // CCSE regtests need special handling of read/write flags
  // for EB checkpoint files
  bool geom_chk_ccse_regtest = false;
  pp.queryAdd("geom_chk_ccse_regtest", geom_chk_ccse_regtest);

  if (geom_chk_ccse_regtest) {

    if (restart_file.empty()) {
      geom_chk_read = false;
      geom_chk_write = true;
    } else {
      geom_chk_read = true;
      geom_chk_write = false;
    }
  }

  if (geom_chk_read) {

    // Override any geometry defined.
    geom_type = "restart";

    // Avoid ambiguous inputs
    if (geom_chk_write) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Both mfix.geom_chk_read and mfix.geom_chk_write cannot be"
          << "\nenabled for the same run";
    }
  }

  // Geometry checkpoint file
  std::string geom_chk_file {"geom_chk"};
  pp.query("geom_chk_file", geom_chk_file);

  // Levelset checkpoint file
  std::string geom_levelset_chk_file {"geom_levelset_chk"};
  pp.query("geom_levelset_chk_file", geom_levelset_chk_file);

  if (geom_type == "box") {

    Print() << "\n Building box geometry.\n";
    make_eb_box(a_geom);
    m_contains_ebs = true;

  } else if (geom_type == "cylinder") {

    Print() << "\n Building cylinder geometry.\n";
    make_eb_cylinder(a_geom);
    m_contains_ebs = true;

  } else if (geom_type == "hopper") {

    Print() << "\n Building hopper geometry.\n";
    make_eb_hopper(a_geom);
    m_contains_ebs = true;

  } else if (geom_type == "generic") {

    Print() << "\n Building generic geometry.\n";
    make_eb_generic(a_geom);
    m_contains_ebs = true;

  } else if (geom_type == "csg") {
#ifdef CSG_EB
    Print() << "\n Building geometry from csg file.\n";
    make_eb_csg(a_geom);
    m_contains_ebs = true;
#else
    reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "CSG geometry selected but the solver was NOT built with CSG support!";
#endif

  } else if (geom_type == "stl") {

    Print() << "\n Building geometry from stl file.\n";
    make_eb_stl(a_geom);
    m_contains_ebs = true;

  } else if (geom_type == "restart") {

    std::ifstream chkptfile(geom_chk_file);
    if (chkptfile.fail()) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "No mfix.geom_chk_file found.";
    }

    if (m_levelset_refinement != 1) {
      std::ifstream refined_chkptfile(geom_levelset_chk_file);
      if (refined_chkptfile.fail()) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "mfix.geom_levelset_chk_file found.";
      }
    }

    Print() << "\n Building geometry from chkptfile: " << geom_chk_file << '\n';
    build_levels_from_chkpt_file(a_geom, geom_chk_file, geom_levelset_chk_file);
    m_contains_ebs = true;

  } else if ( geom_type == "none" ) {
    amrex::Print() << "\n No EB geometry declared.\n";

    make_eb_regular(a_geom);

  } else {

    reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Unknown or invalid input:\n  mfix.geometry: " << geom_type;
  }


  if (geom_chk_write) {
    m_eb_levels[0]->write_to_chkpt_file(geom_chk_file,
        amrex::EB2::ExtendDomainFace(), amrex::EB2::max_grid_size);

    if (m_max_level == 0) {
      if (m_levelset_refinement != 1) {
        m_eb_levels[1]->write_to_chkpt_file(geom_levelset_chk_file,
            amrex::EB2::ExtendDomainFace(), amrex::EB2::max_grid_size);
      }
    }
  }

}


/****************************************************************************
 *                                                                          *
 * Fill EB factories as an initial run. Since the particle container might  *
 * not have been created yet, use mfix grids instead.                       *
 *                                                                          *
 ***************************************************************************/
void MFIXEB::
make_factory ( int const a_lev,
               Geometry            const& a_geom,
               BoxArray            const& a_grids,
               DistributionMapping const& a_dmap)
{
  if (m_debug) { Print() << "Making EB factory on level " << a_lev << '\n'; }

  // Always make m_ebfactory even if all regular.
  m_ebfactory[a_lev] =
        std::make_shared<EBFArrayBoxFactory>(*m_eb_levels[a_lev], a_geom,
        a_grids, a_dmap, nghost_factory() , m_support_level);

  // The particle factory will point to the fluid factory for single-grid
  // runs. We only make the particle factory for dual grid runs.
  m_particle_ebfactory[a_lev] = m_ebfactory[a_lev];
}


void MFIXEB::
update_factory (int a_lev, Geometry const& a_geom,
                BoxArray            const& a_grids,
                DistributionMapping const& a_dmap)
{
  if (m_debug) { Print() << "Updating EB factory on level " << a_lev << '\n'; }

  // Verify that some kind of EB geometry has already been defined
  AMREX_ASSERT(!EB2::IndexSpace::empty());

  m_ebfactory[a_lev] = std::make_shared<EBFArrayBoxFactory>(
      *(m_eb_levels[a_lev]), a_geom, a_grids, a_dmap,
      nghost_factory(), m_support_level);

#if 0
  std::unique_ptr<FabFactory<FArrayBox> > new_factory = makeEBFabFactory(
      a_geom, a_grids, a_dmap, nghost_factory(), EBSupport::full);

  m_factory[a_lev] = std::move(new_factory);
#endif
}


void MFIXEB::
make_particle_factory ( int const a_lev,
                        Geometry            const& a_geom,
                        BoxArray            const& a_grids,
                        DistributionMapping const& a_dmap)
{
  // TODO: This should be changed so that we only create the particle factory
  // if the grids and dmap differ from the fluids. I think we would have
  // to change this to a shared pointer.
  m_particle_ebfactory[a_lev] =
      std::make_shared<EBFArrayBoxFactory>(*m_eb_levels[a_lev], a_geom,
      a_grids,a_dmap, nghost_particle_factory(), m_support_level);
}


void MFIXEB::regrid_levelset_array (int a_lev,
                                    Vector<Geometry> const& a_geom,
                                    MFIXParticleContainer const* a_pc)
{
  if (m_debug) { Print() << "Regridding levelset on level " << a_lev << '\n'; }

  // Verify that some kind of EB geometry has already been defined
   AMREX_ASSERT(!EB2::IndexSpace::empty());

   const DistributionMapping& dm = a_pc->ParticleDistributionMap(a_lev);
   const BoxArray&            ba = a_pc->ParticleBoxArray(a_lev);

   int update = false;
   if ( m_particle_ebfactory[a_lev] == nullptr ) { update = true; }
   else {

     const DistributionMapping&  eb_dm = m_particle_ebfactory[a_lev]->DistributionMap();
     const BoxArray&             eb_ba = m_particle_ebfactory[a_lev]->boxArray();

     if ( (dm != eb_dm) || (ba != eb_ba) ) { update = true; }
   }

   if (update) {

     amrex::Print() << "Updating particle ebfactory\n";

     m_particle_ebfactory[a_lev] =
       std::make_shared<EBFArrayBoxFactory>(*(m_eb_levels[a_lev]), a_geom[a_lev],
         ba, dm, nghost_particle_factory(), m_support_level);

     Print() << "Regridding level-set on lev = " << a_lev << std::endl;

     const BoxArray nd_ba = amrex::convert(ba, IntVect::TheNodeVector());

     std::unique_ptr<MultiFab> new_level_set{new MultiFab()};

     if (m_level_sets[a_lev]->boxArray() == nd_ba) {

       MFUtil::regrid(*new_level_set, nd_ba, dm, *(m_level_sets[a_lev]), true);

     } else {

       int nc = m_level_sets[a_lev]->nComp();
       int ng = m_level_sets[a_lev]->nGrow();
       const Periodicity& period = a_geom[a_lev].periodicity();
       new_level_set->define(nd_ba, dm, nc, ng);
       new_level_set->setVal(0.);

       new_level_set->ParallelCopy(*(m_level_sets[a_lev]), 0, 0, nc, ng, ng, period);
     }

     std::swap(m_level_sets[a_lev], new_level_set);

     //________________________________________________________________________
     // If we're operating in single-level mode, the level-set has a second
     // (refined) MultiFab that also needs to be regridded.

     if ( (m_max_level == 0) && (a_lev == 0)) {

       Print() << "Also regridding refined level-set" << std::endl;

       BoxArray ref_nd_ba = amrex::convert(ba, IntVect::TheNodeVector());
       ref_nd_ba.refine(m_levelset_refinement);

       std::unique_ptr<MultiFab> new_level_set_lev {new MultiFab()};

       if (m_level_sets[a_lev+1]->boxArray() == ref_nd_ba)
       {
         MFUtil::regrid(*new_level_set_lev, ref_nd_ba, dm, *(m_level_sets[a_lev+1]), true);
       }
       else
       {
         int nc = m_level_sets[a_lev+1]->nComp();
         int ng = m_level_sets[a_lev+1]->nGrow();
         const Periodicity& period = a_geom[a_lev].periodicity();
         new_level_set_lev->define(ref_nd_ba, dm, nc, ng);
         new_level_set_lev->setVal(0.0);

         new_level_set_lev->ParallelCopy(*(m_level_sets[a_lev+1]), 0, 0, nc, ng, ng, period);
       }

       std::swap(m_level_sets[a_lev+1], new_level_set_lev);
     }
   }
}


void MFIXEB::
fill_levelsets ( Vector<Geometry> const& a_geom,
                 MFIXParticleContainer const* a_pc,
                 MFIXBoundaryConditions const& a_bcs,
                 MFIXPorousMedia const& a_porous_media)
{
  if (m_debug) { Print() << "fill_levelsets\n"; }
  /****************************************************************************
   *                                                                          *
   * Fill levels either as a single level with refinement, or as a proper     *
   * multi-level hierarchy                                                    *
   *                                                                          *
   ***************************************************************************/

  if (m_max_level == 0) {

    //_______________________________________________________________________
    // COMPATIBILITY: in order to be compatible with benchmarks, fill using
    // the level-set factory (which is then thrown away).

    const DistributionMapping & part_dm = a_pc->ParticleDistributionMap(0);
    const BoxArray &            part_ba = a_pc->ParticleBoxArray(0);

    const DistributionMapping& ls_dm = part_dm;
    const Geometry& ls_geom = amrex::refine(a_geom[0], m_levelset_refinement);

    BoxArray ls_ba = amrex::convert(part_ba, IntVect::TheNodeVector());
    m_level_sets[0] = std::make_unique<MultiFab>(ls_ba, ls_dm, 1,
        m_levelset_pad/m_levelset_refinement);

    if (m_levelset_refinement != 1) ls_ba.refine(m_levelset_refinement);
    m_level_sets[1] = std::make_unique<MultiFab>(ls_ba, ls_dm, 1, m_levelset_pad);

    //___________________________________________________________________________
    // NOTE: Boxes are different (since we're not refining, we need to treat
    // corners this way). IMPORTANT: the box case is assembled from planes
    // => fill the level-set with these (otherwise corners will be
    // inaccurately resolved.)

    ParmParse pp("mfix");

    std::string geom_type;
    pp.query("geometry", geom_type);

    if (geom_type == "box") {

      if ( a_geom[0].isAllPeriodic() ) {

        make_eb_regular(a_geom);
        m_level_sets[1]->setVal(std::numeric_limits<Real>::max());
        m_level_sets[0]->setVal(std::numeric_limits<Real>::max());

      } else {

        const auto plo = a_geom[0].ProbLoArray();
        const auto phi = a_geom[0].ProbHiArray();

        Vector<Real> boxLo(3), boxHi(3);
        Real offset    = 1.0e-15;

        for (int idim(0); idim < 3; ++idim) {

            boxLo[idim] = plo[idim];
            boxHi[idim] = phi[idim];
        }

        ParmParse pp_box("box");

        pp_box.queryarr("Lo", boxLo,  0, 3);
        pp_box.queryarr("Hi", boxHi,  0, 3);

        pp_box.queryAdd("offset", offset);

        Real xlo = ((!a_bcs.ls_wall_flag().test(0)) ? boxLo[0] : plo[0]) + offset;
        Real xhi = ((!a_bcs.ls_wall_flag().test(1)) ? boxHi[0] : phi[0]) - offset;

        Real ylo = ((!a_bcs.ls_wall_flag().test(2)) ? boxLo[1] : plo[1]) + offset;
        Real yhi = ((!a_bcs.ls_wall_flag().test(3)) ? boxHi[1] : phi[1]) - offset;

        Real zlo = ((!a_bcs.ls_wall_flag().test(4)) ? boxLo[2] : plo[2]) + offset;
        Real zhi = ((!a_bcs.ls_wall_flag().test(5)) ? boxHi[2] : phi[2]) - offset;

        RealArray arrLo {xlo, ylo, zlo};
        RealArray arrHi {xhi, yhi, zhi};

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        for (int idim(0); idim<3; ++idim) {
          if (a_geom[0].isPeriodic(idim)) {
            arrLo[idim] = 2.0*plo[idim] - phi[idim];
            arrHi[idim] = 2.0*phi[idim] - plo[idim];
          }
        }

        reporter::Log(reporter::Status)
            << "Bounding box for levelset imposed walls (EB box geometry):"
            << "\n   Lo: " << arrLo[0] << "  " << arrLo[1] << "  " << arrLo[2]
            << "\n   Hi: " << arrHi[0] << "  " << arrHi[1] << "  " << arrHi[2];

        auto if_box = EB2::BoxIF( arrLo, arrHi, /*inside=*/true);
        auto gshop = EB2::makeShop(if_box);

        amrex::FillImpFunc(*m_level_sets[1], gshop, ls_geom);
        m_level_sets[1]->negate(m_level_sets[1]->nGrow()); // signed distance f = - imp. f.
        if (m_levelset_refinement == 1) {
          MultiFab::Copy(*m_level_sets[0], *m_level_sets[1], 0, 0, 1, m_level_sets[0]->nGrow());
        } else {
          amrex::average_down_nodal(*m_level_sets[1], *m_level_sets[0],
                                    IntVect(m_levelset_refinement),
                                    m_level_sets[0]->nGrow(), true);
        }
      }

      return;
    }

    if (m_contains_ebs) {
      amrex::FillSignedDistance(*m_level_sets[1], *m_eb_levels[1], *m_particle_ebfactory[0],
                                m_levelset_refinement);
      if (m_levelset_refinement == 1) {
        MultiFab::Copy(*m_level_sets[0], *m_level_sets[1], 0, 0, 1, m_level_sets[0]->nGrow());
      } else {
        amrex::average_down_nodal(*m_level_sets[1], *m_level_sets[0],
                                  IntVect(m_levelset_refinement),
                                  m_level_sets[0]->nGrow(), true);
      }
    } else {
      m_level_sets[1]->setVal(std::numeric_limits<Real>::max());
      m_level_sets[0]->setVal(std::numeric_limits<Real>::max());
    }
  }
  else
  {
    amrex::Abort("xxxxx fill_levelsets todo");
#if 0
    const DistributionMapping & part_dm = a_pc->ParticleDistributionMap(0);
    const BoxArray &            part_ba = a_pc->ParticleBoxArray(0);

    //_______________________________________________________________________
    // Multi-level level-set: build finer level using coarse level set

    EBFArrayBoxFactory eb_factory(*m_eb_levels[0], a_geom[0], part_ba, part_dm,
                                  {m_levelset_eb_pad + 2, m_levelset_eb_pad + 2,
                                   m_levelset_eb_pad + 2}, EBSupport::full);

    // NOTE: reference BoxArray is not nodal
    BoxArray ba = amrex::convert(part_ba, IntVect::TheNodeVector());
    m_level_sets[0] = std::make_unique<MultiFab>(ba, part_dm, 1, m_levelset_pad);
    iMultiFab valid(ba, part_dm, 1, m_levelset_pad);

    MultiFab impfunc(ba, part_dm, 1, m_levelset_pad);
    m_eb_levels[0]->fillLevelSet(impfunc, a_geom[0]);
    impfunc.FillBoundary(a_geom[0].periodicity());


    LSFactory::fill_data(*m_level_sets[0], valid, *m_particle_ebfactory[0], impfunc,
                         32, 1, 1, a_geom[0], a_geom[0]);

    for (int lev = 1; lev <= m_max_level; lev++)
    {
        const DistributionMapping & part_dm_lev = a_pc->ParticleDistributionMap(lev);
        const BoxArray &            part_ba_lev = a_pc->ParticleBoxArray(lev);

        // NOTE: reference BoxArray is not nodal
        BoxArray ba_lev = amrex::convert(part_ba_lev, IntVect::TheNodeVector());
        if (m_level_sets[lev] != nullptr) delete m_level_sets[lev];
        m_level_sets[lev] = new MultiFab();
        // iMultiFab valid_lev(ba_lev, part_dm_lev, 1, levelset_pad);

        // Fills level-set[lev] with coarse data
        LSCoreBase::MakeNewLevelFromCoarse( *m_level_sets[lev], *m_level_sets[lev-1],
                                           part_ba_lev, part_dm_lev, a_geom[lev], a_geom[lev-1],
                                           bcs_ls, refRatio(lev-1));

        EBFArrayBoxFactory eb_factory_lev(*m_eb_levels[lev], a_geom[lev], part_ba_lev, part_dm_lev,
                                          {m_levelset_eb_pad + 2, m_levelset_eb_pad + 2,
                                           m_levelset_eb_pad + 2}, EBSupport::full);

        MultiFab impfunc_lev(ba_lev, part_dm_lev, 1, m_levelset_pad);
        m_eb_levels[lev]->fillLevelSet(impfunc_lev, a_geom[lev]);
        impfunc_lev.FillBoundary(a_geom[lev].periodicity());

        IntVect ebt_size{32, 32, 32}; // Fudge factors...
        LSCoreBase::FillLevelSet(*m_level_sets[lev], *m_level_sets[lev], eb_factory_lev, impfunc_lev,
                                 ebt_size, m_levelset_eb_pad, a_geom[lev]);
    }
#endif
  }

  // Add walls (for instance MI) to levelset data
  intersect_ls_walls(a_geom, a_bcs, a_porous_media);
}

void MFIXEB::
intersect_ls_walls ( Vector<Geometry> const& a_geom,
                     MFIXBoundaryConditions const& a_bcs,
                     MFIXPorousMedia const& a_porous_media) const
{
  // Skip this routine if there are no LS walls to impose
  if ( a_bcs.ls_wall_flag().none() ) { return; }

  const auto plo = a_geom[0].ProbLoArray();
  const auto phi = a_geom[0].ProbHiArray();

  Real xlo = a_bcs.ls_wall_flag().test(0) ? plo[0] + 1.0e-15 : 2.0*plo[0] - phi[0];
  Real xhi = a_bcs.ls_wall_flag().test(1) ? phi[0] - 1.0e-15 : 2.0*phi[0] - plo[0];

  Real ylo = a_bcs.ls_wall_flag().test(2) ? plo[1] + 1.0e-15 : 2.0*plo[1] - phi[1];
  Real yhi = a_bcs.ls_wall_flag().test(3) ? phi[1] - 1.0e-15 : 2.0*phi[1] - plo[1];

  Real zlo = a_bcs.ls_wall_flag().test(4) ? plo[2] + 1.0e-15 : 2.0*plo[2] - phi[2];
  Real zhi = a_bcs.ls_wall_flag().test(5) ? phi[2] - 1.0e-15 : 2.0*phi[2] - plo[2];

  RealArray arrLo {xlo, ylo, zlo};
  RealArray arrHi {xhi, yhi, zhi};

  reporter::Log(reporter::Status)
      << "Bounding box for levelset imposed walls:"
      << "\n   Lo: " << arrLo[0] << "  " << arrLo[1] << "  " << arrLo[2]
      << "\n   Hi: " << arrHi[0] << "  " << arrHi[1] << "  " << arrHi[2];

  auto bounding_box = EB2::BoxIF( arrLo, arrHi, /*inside=*/true);

  auto gshop = EB2::makeShop(bounding_box);

  if (m_max_level == 0)
  {
    //_______________________________________________________________________
    // Baseline Level-Set
    {
      const int ng = m_level_sets[0]->nGrow();
      const BoxArray & ba = m_level_sets[0]->boxArray();
      const DistributionMapping & dm = m_level_sets[0]->DistributionMap();
      MultiFab wall_if(ba, dm, 1, ng);

      amrex::FillImpFunc(wall_if, gshop, a_geom[0]);

      // The difference between this and the previous code is we no
      // longer clamp the level sets.
      for (MFIter mfi(*m_level_sets[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox();
        Array4<Real> const& sdf = m_level_sets[0]->array(mfi);
        Array4<Real const> const& wif = wall_if.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          sdf(i,j,k) = amrex::min(sdf(i,j,k),-wif(i,j,k));
        });
      }

      if (a_porous_media.nregions() > 0) {
        a_porous_media.block_particles(*m_level_sets[0], a_geom[0]);
      }
    }

    //_______________________________________________________________________
    // Refined Level-Set
    if (m_levelset_refinement == 1) {
      MultiFab::Copy(*m_level_sets[1], *m_level_sets[0], 0, 0, 1, m_level_sets[0]->nGrow());
    } else {
      const int ng = m_level_sets[1]->nGrow();
      const BoxArray & ba = m_level_sets[1]->boxArray();
      const DistributionMapping & dm = m_level_sets[1]->DistributionMap();
      MultiFab wall_if(ba, dm, 1, ng);

      amrex::FillImpFunc(wall_if, gshop, amrex::refine(a_geom[0], m_levelset_refinement));

      // The difference between this and the previous code is we no
      // longer clamp the level sets.
      for (MFIter mfi(*m_level_sets[1],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox();
        Array4<Real> const& sdf = m_level_sets[1]->array(mfi);
        Array4<Real const> const& wif = wall_if.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          sdf(i,j,k) = amrex::min(sdf(i,j,k),-wif(i,j,k));
        });
      }

      if (a_porous_media.nregions() > 0) {
        a_porous_media.block_particles(*m_level_sets[1],
            amrex::refine(a_geom[0], m_levelset_refinement));
      }
    }
  }
  else
  {
    amrex::Abort("xxxxx intersect_ls_walls todo");
#if 0
    //_______________________________________________________________________
    // Multi-level level-set: apply wall-intersection to each level

    for (int lev = 0; lev <= m_max_level; lev++)
    {
      const int ng = m_level_sets[lev]->nGrow();
      const BoxArray & ba = m_level_sets[lev]->boxArray();
      const DistributionMapping & dm = m_level_sets[lev]->DistributionMap();

      MultiFab wall_if(ba, dm, 1, ng);
      iMultiFab valid_lev(ba, dm, 1, ng);
      valid_lev.setVal(1);

      GShopLSFactory<UnionListIF<EB2::PlaneIF>> gshop_lsf(gshop, a_geom[lev], ba, dm, ng);
      std::unique_ptr<MultiFab> impfunc = gshop_lsf.fill_impfunc();

      LSFactory::fill_data(wall_if, valid_lev, *impfunc, m_levelset_eb_pad, a_geom[lev]);
      LSFactory::intersect_data(*m_level_sets[lev], valid_lev, wall_if, valid_lev, a_geom[lev]);

      if (a_porous_media.nregions() > 0) {
        a_porous_media.block_particles(*m_level_sets[lev], a_geom[lev]);
      }
    }
#endif
  }
}

void MFIXEB::
build_levels_from_chkpt_file ( Vector<Geometry> const& a_geom,
                               std::string const& geom_chk_file,
                               std::string const& geom_levelset_chk_file) {

  EB2::BuildFromChkptFile(geom_chk_file, a_geom[m_max_level], m_max_level, 100);

  const EB2::IndexSpace& ebis = EB2::IndexSpace::top();
  for (int lev = 0; lev <= m_max_level; lev ++) {

    m_eb_levels[lev] = &(ebis.getLevel(a_geom[lev]));
  }

  if (m_max_level == 0) {

    if (m_levelset_refinement == 1) {

      m_eb_levels[1] = m_eb_levels[0];

    } else {

      Geometry geom_ls = amrex::refine(a_geom[0], m_levelset_refinement);
      EB2::BuildFromChkptFile(geom_levelset_chk_file, geom_ls, 0, 100);
      m_eb_levels[1] = &(EB2::IndexSpace::top().getLevel(geom_ls));
    }
  }
}
