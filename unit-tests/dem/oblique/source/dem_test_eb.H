#ifndef DEM_TEST_EB_H_
#define DEM_TEST_EB_H_

#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>

#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>

class dem_test_eb {

  public:

    dem_test_eb ( amrex::Vector<amrex::Geometry> const& a_geom, int const a_max_lev);

    void make_factory ( amrex::Vector<amrex::Geometry>            const& a_geom,
                        amrex::Vector<amrex::DistributionMapping> const& a_dmap,
                        amrex::Vector<amrex::BoxArray>            const& a_grids);

    int nghost_basic () const  { return 2; }
    int nghost_volume () const { return 2; }
    int nghost_full () const   { return 2; }

    amrex::Vector< std::unique_ptr<amrex::EBFArrayBoxFactory>> m_factory;

    amrex::Vector< const amrex::FabArray<amrex::EBCellFlagFab>* > flags () noexcept
    { amrex::Vector< const amrex::FabArray<amrex::EBCellFlagFab>* > r;
      r.reserve(m_factory.size());
      for (auto& lev_factory : m_factory) {
        r.push_back( &(lev_factory->getMultiEBCellFlagFab()) );
      }
      return r;
    }

    amrex::Vector< const amrex::MultiFab* > volfrac () noexcept
    { amrex::Vector< const amrex::MultiFab* > r;
      r.reserve(m_factory.size());
      for (auto& lev_factory : m_factory) {
        r.push_back( &(lev_factory->getVolFrac()) );
      }
      return r;
    }

    amrex::EBFArrayBoxFactory* factory (int a_lev) noexcept
    { return m_factory[a_lev].get(); }

    amrex::Vector< amrex::EBFArrayBoxFactory const*> const_factory () const noexcept
    { return GetVecOfConstPtrs(m_factory); }

  private:

    std::string m_type;

    amrex::EBSupport m_support_level;

    //! EB level constructed from building GeometryShop
    amrex::Vector<const amrex::EB2::Level*> m_eb_levels;

    void make_box (amrex::Vector<amrex::Geometry> const& geom,
                   int const max_level);

    //! Construct EB levels from Geometry shop.
    template<class F>
    void build_level ( amrex::EB2::GeometryShop<F> a_gshop,
                       amrex::Vector<amrex::Geometry> const& a_geom,
                       int const a_max_lev) {

      int const req_lev(100);

      amrex::EB2::Build(a_gshop, a_geom[a_max_lev], a_max_lev, req_lev);

      m_eb_levels.resize(a_max_lev+1);

      const amrex::EB2::IndexSpace& ebis = amrex::EB2::IndexSpace::top();
      for (int lev(0); lev<a_max_lev+1; ++lev) {
        m_eb_levels[lev] = &(ebis.getLevel(a_geom[lev]));
      }
    }

};
#endif
