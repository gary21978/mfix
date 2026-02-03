#ifndef MFIX_POROUS_MEDIA_H_
#define MFIX_POROUS_MEDIA_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>

#include <mfix_regions.H>
#include <mfix_fluid.H>
#include <mfix_run_on.H>

class PorousMediaRegions
{
  public:

    AMREX_GPU_HOST_DEVICE
    PorousMediaRegions ();

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool contains (int lev, int i, int j, int k, int n_pm) const
    {
      return m_box[m_nregions*lev + n_pm].contains(i,j,k);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const amrex::RealBox& realbox (int n_pm) const
    {
      return m_realbox[n_pm];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real vfrac (int n_pm) const
    {
      return m_vfrac[n_pm];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real c1 (int n_pm) const
    {
      return m_c1[n_pm];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real c2 (int n_pm) const
    {
      return m_c2[n_pm];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int allow_particles (int n_pm) const
    {
      return m_allow_particles[n_pm];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool contains (int lev, int i, int j, int k) const
    {
      for (int n = 0; n < m_nregions; ++n) {
        if (contains(lev,i,j,k,n)) { return true; }
      }

      return false;
    }

  private:
    friend class MFIXPorousMedia;

    int m_nregions;
    amrex::Box* m_box;
    amrex::RealBox* m_realbox;
    amrex::Real* m_vfrac;
    amrex::Real* m_c1;
    amrex::Real* m_c2;
    int* m_allow_particles;
};

class MFIXPorousMedia
{
  public:
    // constructor
    MFIXPorousMedia ();

    // destructor
    ~MFIXPorousMedia ();

    void Initialize (const MFIXRegions& regions,
                     const amrex::Vector<amrex::Geometry>& geom);

    int nregions () const;

    void impose_vfrac (amrex::Vector<amrex::MultiFab*> const& ep_g) const;

    void ComputeVelSource ( amrex::Vector< amrex::MultiFab      *> const& a_S_p,
                            amrex::Vector< amrex::MultiFab      *> const& a_S_c,
                            amrex::Vector< amrex::MultiFab const*> const& a_rho,
                            amrex::Vector< amrex::MultiFab const*> const& a_vel,
                            amrex::Vector< amrex::MultiFab const*> const& a_T,
                            amrex::Vector< amrex::MultiFab const*> const& a_X,
                            ThermoPropertyData const& fluid_props,
                            amrex::Real a_S_fac) const;

    void block_particles (amrex::MultiFab& a_level_set,
                          amrex::Geometry const& a_geom) const;

    bool intersects_box (int lev, amrex::Box const& box) const;

    template <amrex::RunOn run_on>
    const PorousMediaRegions& regions () const
    {
      if constexpr (run_on == amrex::RunOn::Cpu) {
        AMREX_ALWAYS_ASSERT(m_h_pm_regions);
        return *m_h_pm_regions;
      } else if constexpr (run_on == amrex::RunOn::Gpu) {
        AMREX_ALWAYS_ASSERT(m_d_pm_regions);
        return *m_d_pm_regions;
      }
    }

  private:
    int m_nregions;

    PorousMediaRegions* m_h_pm_regions;
    PorousMediaRegions* m_d_pm_regions;

    amrex::Gpu::HostVector<amrex::Box> m_box;
    amrex::Gpu::DeviceVector<amrex::Box> m_d_box;

    amrex::Gpu::HostVector<amrex::RealBox> m_realbox;
    amrex::Gpu::DeviceVector<amrex::RealBox> m_d_realbox;

    amrex::Gpu::HostVector<amrex::Real> m_vfrac;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_vfrac;

    amrex::Gpu::HostVector<amrex::Real> m_c1;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_c1;

    amrex::Gpu::HostVector<amrex::Real> m_c2;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_c2;

    amrex::Gpu::HostVector<int> m_allow_particles;
    amrex::Gpu::DeviceVector<int> m_d_allow_particles;
};

#endif
