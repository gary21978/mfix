#include <post_mfix.H>
#include <pm_diffusion.H>

using namespace amrex;

Real
post_mfix::
calc_fluct ( MultiFab const* a_MF,
             int const a_comp, Real const a_avg )
{

  // Step 2: Compute fluct = sum((ep_g(i,j,k) - <ep_g>)^2)
  // ---------------------------------------------------------------------------

  auto const& xma = a_MF->const_arrays();
  Real sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, *a_MF, IntVect(0),
  [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
  {
    auto const& xfab = xma[box_no];
    Real t = xfab(i,j,k,a_comp) - a_avg;
    return t*t;
  });

  return sm;
}
