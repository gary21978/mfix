#ifndef MFIX_UTIL_H
#define MFIX_UTIL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Box.H>
#include <AMReX_Vector.H>
#include <AMReX_EBFabFactory.H>


struct Utils {


static amrex::IntVect
locate_max_eps (amrex::Vector< amrex::MultiFab const*> const& ep_s_in,
                amrex::Real const max_eps)
{

  int imax(-100), jmax(-100), kmax(-100);

#ifdef AMREX_USE_GPU
  amrex::Gpu::DeviceScalar<int> imax_gpu(-100);
  amrex::Gpu::DeviceScalar<int> jmax_gpu(-100);
  amrex::Gpu::DeviceScalar<int> kmax_gpu(-100);

  int* p_imax = imax_gpu.dataPtr();
  int* p_jmax = jmax_gpu.dataPtr();
  int* p_kmax = kmax_gpu.dataPtr();
#endif

  for (int lev(0); lev < ep_s_in.size(); lev++) {

    constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

    for (amrex::MFIter mfi(*ep_s_in[lev],false); mfi.isValid(); ++mfi) {
      amrex::Box const& bx = mfi.tilebox();
      // Array4<Real const> const& epg = ld.ep_g->const_array(mfi);
      amrex::Array4<const amrex::Real> const& eps = ep_s_in[lev]->const_array(mfi);

#ifdef AMREX_USE_GPU
      amrex::ParallelFor(bx, [eps, max_eps, p_imax, p_jmax, p_kmax]
#else
      amrex::ParallelFor(bx, [eps, max_eps, &imax, &jmax, &kmax]
#endif
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if( amrex::Math::abs(eps(i,j,k) - max_eps) < tolerance ){
#ifdef AMREX_USE_GPU
          *p_imax = i;
          *p_jmax = j;
          *p_kmax = k;
#else
          imax = i;
          jmax = j;
          kmax = k;
#endif
        }
      });
    } // mfi
  } // lev


  amrex::Gpu::synchronize();

#ifdef AMREX_USE_GPU
  imax = imax_gpu.dataValue();
  jmax = jmax_gpu.dataValue();
  kmax = kmax_gpu.dataValue();
#endif

  amrex::ParallelDescriptor::ReduceIntMax({imax, jmax, kmax});

  return amrex::IntVect{imax, jmax, kmax};
}


};

#endif
