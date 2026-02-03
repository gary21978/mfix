#include <AMReX_ParmParse.H>
#include <AMReX_EB_utils.H>
#include <AMReX_EB2_IF_Box.H>

#include <mfix_porous_media.H>
#include <mfix_reporter.H>
#include <mfix_run_on.H>

#include <cmath>

using namespace amrex;

AMREX_GPU_HOST_DEVICE
PorousMediaRegions::PorousMediaRegions ()
  : m_nregions(0)
  , m_box(nullptr)
  , m_realbox(nullptr)
  , m_vfrac(nullptr)
  , m_c1(nullptr)
  , m_c2(nullptr)
  , m_allow_particles(nullptr)
{}

MFIXPorousMedia::MFIXPorousMedia ()
  : m_nregions(0)
  , m_h_pm_regions(nullptr)
  , m_d_pm_regions(nullptr)
{}

MFIXPorousMedia::~MFIXPorousMedia ()
{
  if (m_h_pm_regions != nullptr) {
    delete m_h_pm_regions;
  }

  if (m_d_pm_regions != nullptr) {
    delete m_d_pm_regions;
  }
}

void MFIXPorousMedia::
Initialize ( const MFIXRegions& regions,
             const Vector<Geometry>& geom)
{
  const int nlev = geom.size();

  ParmParse pp("pm");

  std::vector<std::string> input_regions;
  pp.queryarr("regions", input_regions);

  // Nothing to do if there are no porous media regions
  if (input_regions.size() == 0) {
    return;
  }

  // Loop over porous media inputs
  for (size_t i = 0; i < input_regions.size(); ++i) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(regions.getIndex(input_regions[i]) != -1,
          "Invalid pm region");
    ++m_nregions;
  }

  AMREX_ALWAYS_ASSERT(static_cast<int>(input_regions.size()) == m_nregions);

  m_box.resize(m_nregions*nlev);
  m_realbox.resize(m_nregions);
  m_vfrac.resize(m_nregions);
  m_c1.resize(m_nregions);
  m_c2.resize(m_nregions);
  m_allow_particles.resize(m_nregions);

  for (size_t i = 0; i < input_regions.size(); ++i) {
    std::string pm_region_prefix = "pm."+input_regions[i];
    amrex::ParmParse ppPM(pm_region_prefix);

    Real lvfrac = 0., lc1 = 0., lc2 = 0.;
    int lallow = 0;

    if (!ppPM.query("volfrac", lvfrac)) {
       reporter::Log(reporter::Error,__FILE__, __LINE__)
         << "Failed to process porous media volfrac.\n"
         << "PM region: " << input_regions[i] << '\n'
         << "Please correct the input deck.";
    }

    if (!ppPM.query("c1", lc1)) {
       reporter::Log(reporter::Error,__FILE__, __LINE__)
         << "Failed to process porous media c1.\n"
         << "PM region: " << input_regions[i] << '\n'
         << "Please correct the input deck.";
    }

    if (!ppPM.query("c2", lc2)) {
       reporter::Log(reporter::Error,__FILE__, __LINE__)
         << "Failed to process porous media c2.\n"
         << "PM region: " << input_regions[i] << '\n'
         << "Please correct the input deck.";
    }

    ppPM.query("allow_particles", lallow);

    // Compute box from realbox
    // Also check for grid-alignment
    const auto& rb = *regions.getRegion(input_regions[i]);

    for (int lev = 0; lev < nlev; ++lev) {
      const auto& dx = geom[lev].CellSize();
      const auto& plo = geom[lev].ProbLo();
      IntVect small, big;

      for (int dir = 0; dir < 3; ++dir) {
        Real rlo = (rb.lo(dir)-plo[dir])/dx[dir];
        Real rhi = (rb.hi(dir)-plo[dir])/dx[dir];
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
               almostEqual(rlo, std::round(rlo), 10)
            && almostEqual(rhi, std::round(rhi), 10),
            "Porous media region must be grid-aligned");

        small[dir] = static_cast<int>(std::round(rlo));
        big[dir] = static_cast<int>(std::round(rhi)) - 1;
      }

      m_box[m_nregions*lev + i] = Box(small, big);
    }

    m_realbox[i] = rb;
    m_vfrac[i] = lvfrac;
    m_c1[i] = lc1;
    m_c2[i] = lc2;
    m_allow_particles[i] = lallow;
  }

  // Allocate m_h regions
  m_h_pm_regions = new PorousMediaRegions();

  m_h_pm_regions->m_nregions = m_nregions;

  if (m_nregions > 0) {
    m_h_pm_regions->m_box = m_box.data();
    m_h_pm_regions->m_realbox = m_realbox.data();
    m_h_pm_regions->m_vfrac = m_vfrac.data();
    m_h_pm_regions->m_c1 = m_c1.data();
    m_h_pm_regions->m_c2 = m_c2.data();
    m_h_pm_regions->m_allow_particles = m_allow_particles.data();
  }

  // Allocate m_d regions
#ifdef AMREX_USE_GPU
  m_d_pm_regions = new PorousMediaRegions();

  m_d_pm_regions->m_nregions = m_nregions;

  if (m_nregions > 0) {
    m_d_box.resize(m_nregions*nlev);
    m_d_realbox.resize(m_nregions);
    m_d_vfrac.resize(m_nregions);
    m_d_c1.resize(m_nregions);
    m_d_c2.resize(m_nregions);
    m_d_allow_particles.resize(m_nregions);

    Gpu::copyAsync(Gpu::hostToDevice, m_box.begin(), m_box.end(),
          m_d_box.begin());
    m_d_pm_regions->m_box = m_d_box.data();

    Gpu::copyAsync(Gpu::hostToDevice, m_realbox.begin(), m_realbox.end(),
          m_d_realbox.begin());
    m_d_pm_regions->m_realbox = m_d_realbox.data();

    Gpu::copyAsync(Gpu::hostToDevice, m_vfrac.begin(), m_vfrac.end(),
          m_d_vfrac.begin());
    m_d_pm_regions->m_vfrac = m_d_vfrac.data();

    Gpu::copyAsync(Gpu::hostToDevice, m_c1.begin(), m_c1.end(),
          m_d_c1.begin());
    m_d_pm_regions->m_c1 = m_d_c1.data();

    Gpu::copyAsync(Gpu::hostToDevice, m_c2.begin(), m_c2.end(),
          m_d_c2.begin());
    m_d_pm_regions->m_c2 = m_d_c2.data();

    Gpu::copyAsync(Gpu::hostToDevice, m_allow_particles.begin(), m_allow_particles.end(),
          m_d_allow_particles.begin());
    m_d_pm_regions->m_allow_particles = m_d_allow_particles.data();
  }
#endif

}

int MFIXPorousMedia::nregions () const
{
  return m_nregions;
}

void MFIXPorousMedia::impose_vfrac (Vector<MultiFab*> const& ep_g) const
{
  const int nlev = ep_g.size();

  int n_pm_regions = m_nregions;
#ifdef AMREX_USE_GPU
  const auto& pm_regions = *m_d_pm_regions;
#else
  const auto& pm_regions = *m_h_pm_regions;
#endif

  for (int lev = 0; lev < nlev; lev++) {
    for (MFIter mfi(*ep_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      // Tilebox
      Box bx = mfi.tilebox();
      Array4<Real> const& ep_g_arr  = ep_g[lev]->array(mfi);

      ParallelFor(bx,[ep_g_arr,n_pm_regions,pm_regions,lev]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        for (int n = 0; n < n_pm_regions; ++n) {
          if (pm_regions.contains(lev,i,j,k,n)) {
            if (pm_regions.allow_particles(n) == 0) {
              ep_g_arr(i,j,k) = 1. - pm_regions.vfrac(n);
            } else {
              ep_g_arr(i,j,k) -= pm_regions.vfrac(n);
            }
          }
        }
      });
    }
  }
}

void MFIXPorousMedia::
ComputeVelSource ( Vector< MultiFab      *> const& a_S_p,
                   Vector< MultiFab      *> const& a_S_c,
                   Vector< MultiFab const*> const& a_rho,
                   Vector< MultiFab const*> const& a_vel,
                   Vector< MultiFab const*> const& a_T,
                   Vector< MultiFab const*> const& a_X,
                   ThermoPropertyData const& fluid_props,
                   Real a_S_fac) const
{
  int n_pm_regions = m_nregions;

  if (n_pm_regions < 1) { return; }

#ifdef AMREX_USE_GPU
  const auto& pm_regions = *m_d_pm_regions;
#else
  const auto& pm_regions = *m_h_pm_regions;
#endif

  int const nlev( a_rho.size() );

  Real const S_p_fac(a_S_fac);
  Real const S_c_fac(a_S_fac-1.0);

  for (int lev(0); lev < nlev; ++lev) {

    for (MFIter mfi(*a_rho[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      Box const& bx = mfi.tilebox();

      Array4<Real      > const& S_p = a_S_p[lev]->array(mfi);
      Array4<Real      > const& S_c = a_S_c[lev]->array(mfi);

      Array4<Real const> const& rho = a_rho[lev]->const_array(mfi);
      Array4<Real const> const& vel = a_vel[lev]->const_array(mfi);

      Array4<Real const> const& T = a_T[lev]->const_array(mfi);
      Array4<Real const> const& X = a_X[lev]->const_array(mfi);

      ParallelFor(bx, n_pm_regions, [fluid_props, vel, rho, T, X,
      pm_regions, S_p, S_c, S_p_fac, S_c_fac, lev]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        if (pm_regions.contains(lev,i,j,k,n)) {
          amrex::Real vmag = std::sqrt(vel(i,j,k,0)*vel(i,j,k,0) +
                                       vel(i,j,k,1)*vel(i,j,k,1) +
                                       vel(i,j,k,2)*vel(i,j,k,2));

          Real term = (fluid_props.molViscosity(i,j,k,T,X) / pm_regions.c1(n))
                     +(0.5*pm_regions.c2(n)*rho(i,j,k)*vmag);

          S_p(i,j,k) +=  S_p_fac * term;

          S_c(i,j,k,0) +=  S_c_fac * vel(i,j,k,0) * term;
          S_c(i,j,k,1) +=  S_c_fac * vel(i,j,k,1) * term;
          S_c(i,j,k,2) +=  S_c_fac * vel(i,j,k,2) * term;
        }
      });
    } // MFIter
  } // lev
}

void MFIXPorousMedia::block_particles (MultiFab& a_level_set,
                                       Geometry const& a_geom) const
{
  for (int n = 0; n < m_nregions; ++n) {
    if (m_allow_particles[n] == 0) {
      Real offset = 1.0e-15;
      const auto& rb = m_realbox[n];

      RealArray lo{rb.lo()[0]-offset, rb.lo()[1]-offset, rb.lo()[2]-offset};
      RealArray hi{rb.hi()[0]+offset, rb.hi()[1]+offset, rb.hi()[2]+offset};

      EB2::BoxIF box_if(lo, hi, false);

      auto gshop = EB2::makeShop(box_if);

      const int ng = a_level_set.nGrow();
      const BoxArray & ba = a_level_set.boxArray();
      const DistributionMapping & dm = a_level_set.DistributionMap();
      MultiFab wall_if(ba, dm, 1, ng);

      FillImpFunc(wall_if, gshop, a_geom);

      for (MFIter mfi(a_level_set,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.growntilebox();

        Array4<Real> const& sdf = a_level_set.array(mfi);
        Array4<Real const> const& wif = wall_if.const_array(mfi);

        ParallelFor(bx,[sdf,wif]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          sdf(i,j,k) = amrex::min(sdf(i,j,k),-wif(i,j,k));
        });
      }
    }
  }
}

bool MFIXPorousMedia::intersects_box (int lev, Box const& bx) const
{
  if (m_nregions <= 0) {
    return false;
  }

  for (int n = 0; n < m_nregions; ++n) {
    if (m_box[m_nregions*lev + n].intersects(bx)) {
      return true;
    }
  }

  return false;
}
