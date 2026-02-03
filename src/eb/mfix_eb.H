#ifndef MFIX_EB_H_
#define MFIX_EB_H_

#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>

#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB2_IndexSpace_STL.H>

#include <mfix_pc.H>
#include <mfix_bc.H>
#include <mfix_eb_if.H>
#include <mfix_porous_media.H>

#include <string>

class MFIXEB {

  public:

    MFIXEB (bool debug, int max_level);

    void make_geometry (amrex::Vector<amrex::Geometry> const& a_geom,
                        std::string const& restart_file);

    void make_factory ( int const a_lev,
                        amrex::Geometry            const& a_geom,
                        amrex::BoxArray            const& a_grids,
                        amrex::DistributionMapping const& a_dmap);

    void clear_factory ( int const a_lev )
    { m_ebfactory[a_lev].reset(); }

    //! Remake the EB factory on a_lev
    void update_factory (int a_lev, amrex::Geometry const& a_geom,
                         amrex::BoxArray            const& a_ba,
                         amrex::DistributionMapping const& a_dm);

    void make_particle_factory ( int const a_lev,
                                 amrex::Geometry            const& a_geom,
                                 amrex::BoxArray            const& a_grids,
                                 amrex::DistributionMapping const& a_dmap);


    //! This function regrids the level set function and updates
    //! `particle_ebfactory`. This has to be done separately from the regridding of
    //! the other field variables since LS and `particle_ebfactory` "live" on the
    //! particle grids. Furthermore, the LS always has an additional "refined" level
    //! if not operating in multi-level mode. Hence slightly different regridding
    //! rules are needed.
    void regrid_levelset_array (int a_lev,
                                amrex::Vector<amrex::Geometry> const& a_geom,
                                MFIXParticleContainer const* a_pc);

    // Four ghost cells are required when using EB
    // const int nghost = 4;
    // const int m_eb_basic_grow_cells = nghost;
    // const int m_eb_volume_grow_cells = nghost;
    // const int m_eb_full_grow_cells = nghost;

    int nghost_basic () const  { return 5; }

    int nghost_volume () const { return 5; }

    int nghost_full () const   { return 4; }

    void fill_levelsets (amrex::Vector<amrex::Geometry> const& a_geom,
                         MFIXParticleContainer const* a_pc,
                         MFIXBoundaryConditions const& a_bcs,
                         MFIXPorousMedia const& a_porous_media);

    amrex::Vector<int> nghost_factory () const
    { return {nghost_basic(), nghost_volume(), nghost_full()}; }

    amrex::Vector<int> nghost_particle_factory () const
    { return {amrex::max(nghost_basic(),  m_levelset_pad + 3),
              amrex::max(nghost_volume(), m_levelset_pad + 3),
              amrex::max(nghost_full(),   m_levelset_pad + 3)}; }

    amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory> > const& factory () noexcept
    { return m_ebfactory; }

    amrex::Vector< amrex::EBFArrayBoxFactory const*> factory_const () const noexcept {
      amrex::Vector<amrex::EBFArrayBoxFactory const*> r;
      int nlev(m_ebfactory.size());
      r.reserve( nlev );
      for (int lev(0); lev < nlev; ++lev) {
        r.push_back( m_ebfactory[lev].get() );
      }
      return r;
    }

    amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory> > const& particle_factory () noexcept
    { return m_particle_ebfactory; }


    amrex::Vector< amrex::EBFArrayBoxFactory const*> const pc_factory_const () noexcept {
      amrex::Vector<amrex::EBFArrayBoxFactory const*> r;
      int nlev(m_particle_ebfactory.size());
      r.reserve( nlev );
      for (int lev(0); lev < nlev; ++lev) {
        r.push_back( m_particle_ebfactory[lev].get() );
      }
      return r;
    }

    amrex::EBFArrayBoxFactory const&
    Factory (int lev) const noexcept
    {
      return static_cast<amrex::EBFArrayBoxFactory const&>(*m_ebfactory[lev]);
    }

    amrex::EBFArrayBoxFactory const&
    ParticleFactory (int lev) const noexcept
    {
      return static_cast<amrex::EBFArrayBoxFactory const&>(*m_particle_ebfactory[lev]);
    }

    int levelset_refinement () const noexcept
    { return m_levelset_refinement; }

    int levelset_eb_refinement () const noexcept
    { return m_levelset_eb_refinement; }

    int levelset_pad () const noexcept
    { return m_levelset_pad; }

    int levelset_eb_pad () const noexcept
    { return m_levelset_eb_pad; }

    amrex::Vector< std::unique_ptr<amrex::MultiFab> > const& level_sets () noexcept
    { return m_level_sets; }

    amrex::Vector<const amrex::EB2::Level*> const& levels () noexcept
    { return m_eb_levels; }

    amrex::EBSupport m_support_level;

    void write_surface (amrex::Vector<amrex::Geometry> const& a_geom,
                        amrex::Vector<amrex::DistributionMapping> const& a_dmap,
                        amrex::Vector<amrex::BoxArray> const& a_grids ) const;

    // Public for CUDA
    void intersect_ls_walls ( amrex::Vector<amrex::Geometry> const& a_geom,
                              MFIXBoundaryConditions const& a_bcs,
                              MFIXPorousMedia const& a_porous_media) const;

    void make_eb_box (amrex::Vector<amrex::Geometry> const& a_geom);
    void make_eb_cylinder (amrex::Vector<amrex::Geometry> const& a_geom);
    void make_eb_hopper (amrex::Vector<amrex::Geometry> const& a_geom);
    void make_eb_regular (amrex::Vector<amrex::Geometry> const& a_geom);
    void make_eb_generic (amrex::Vector<amrex::Geometry> const& a_geom);

  private:

    bool m_debug;
    int m_max_level;
    bool m_contains_ebs = false;
    bool m_write_surface = false;

#ifdef CSG_EB
    void make_eb_csg(amrex::Vector<amrex::Geometry> const& a_geom);
#endif
    void make_eb_stl(amrex::Vector<amrex::Geometry> const& a_geom);

    /***************************************************************************
     *                                                                         *
     * EB Data (including level-set data)                                      *
     *                                                                         *
     **************************************************************************/

    //! EB levels representing fluid boundary conditions
    amrex::Vector<const amrex::EB2::Level*> m_eb_levels;

    //! EB factory that lives on the fluid grids
    amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory> > m_ebfactory;

    //! EB factory that lives on the particle grids
    amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory> > m_particle_ebfactory;

    // Level-Set Data => used for particle-wall collisions and fluid
    // reconstruction in particle drag calculation. NOTE: this has at least 2
    // levels: 0) fluid reconstruction, 1..N) particle-wall (refined) collisions
    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_level_sets;

    ///
    /// Parameters describing the level-set grid:
    ///
    int m_levelset_refinement = 1;    // refinement (wrt particle grid) of the level-set's grid
    int m_levelset_eb_refinement = 1; // refinement of the EB facets used to generate the level-set
    int m_levelset_pad = 2;           // padding of the level-set grid
    int m_levelset_eb_pad = 2;        // padding of the EBIS used to generate the level-set (useful for
                                      // ensuring that neighbor grids are considered when filling LS)

    //! Construct EB levels from Geometry shop. This builds each EB level and
    //! saves a pointer to each level into `eb_levels`.
    template<class F> void build_levels (amrex::Vector<amrex::Geometry> const& a_geom,
                                         amrex::EB2::GeometryShop<F> a_gshop)
    {
      /************************************************************************
       *                                                                      *
       * Build EB levels                                                      *
       *                                                                      *
       ***********************************************************************/
      amrex::EB2::Build(a_gshop, a_geom[m_max_level], m_max_level, 100);
      const amrex::EB2::IndexSpace& ebis = amrex::EB2::IndexSpace::top();
      for (int lev = 0; lev <= m_max_level; lev ++) {
        m_eb_levels[lev] = &(ebis.getLevel(a_geom[lev]));
      }

      if (m_max_level == 0)
      {
        if (m_levelset_refinement == 1) {
          m_eb_levels[1] = m_eb_levels[0];
        } else {
          amrex::Geometry geom_ls = amrex::refine(a_geom[0],m_levelset_refinement);
          amrex::EB2::Build(a_gshop, geom_ls, 0, 100);
          m_eb_levels[1] = &(amrex::EB2::IndexSpace::top().getLevel(geom_ls));
        }
      }
    }

    void build_levels (amrex::Vector<amrex::Geometry> const& a_geom,
                       const std::string& stl_file, amrex::Real stl_scale,
                       amrex::Array<amrex::Real,3> const& stl_center, int stl_reverse_normal,
                       bool use_bvh)
    {
      /************************************************************************
       *                                                                      *
       * Build EB levels from STL                                             *
       *                                                                      *
       ***********************************************************************/
      amrex::EB2::IndexSpace::push(new amrex::EB2::IndexSpaceSTL(stl_file, stl_scale,
                                   {stl_center[0], stl_center[1], stl_center[2]},
                                   stl_reverse_normal,
                                   a_geom[m_max_level], m_max_level, 100,
                                   4 /*default ngrow*/, true /*default build_coarse_level_by_coarsening*/,
                                   amrex::EB2::ExtendDomainFace(), amrex::EB2::NumCoarsenOpt(),
                                   use_bvh));

      const amrex::EB2::IndexSpace& ebis = amrex::EB2::IndexSpace::top();
      for (int lev = 0; lev <= m_max_level; lev ++) {
        m_eb_levels[lev] = &(ebis.getLevel(a_geom[lev]));
      }

      if (m_max_level == 0) {

        if (m_levelset_refinement == 1) {

          m_eb_levels[1] = m_eb_levels[0];

        } else {
          amrex::Geometry geom_ls = amrex::refine(a_geom[0], m_levelset_refinement);
          amrex::EB2::IndexSpace::push(new amrex::EB2::IndexSpaceSTL(stl_file, stl_scale,
                                       {stl_center[0], stl_center[1], stl_center[2]},
                                       stl_reverse_normal,
                                       geom_ls, 0, 100,
                                       4 /*default ngrow*/, true /*default build_coarse_level_by_coarsening*/,
                                       amrex::EB2::ExtendDomainFace(), amrex::EB2::NumCoarsenOpt(),
                                       use_bvh));
          m_eb_levels[1] = &(amrex::EB2::IndexSpace::top().getLevel(geom_ls));
        }
      }
    }

    void build_levels_from_chkpt_file (amrex::Vector<amrex::Geometry> const& a_geom,
                                       std::string const& geom_chk_file,
                                       std::string const& geom_levelset_chk_file);

};

#endif
