#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

#include <AMReX_EB_Redistribution.H>

#include <mfix_deposition_op.H>
#include <mfix_mf_helpers.H>

using namespace amrex;

void MFIXDepositionOp::
Redistribute ( int const /*a_lev*/)
{
  int const lev0(0);

  if (m_redist_type == m_MaxPack) {

    Redistribute_MaxPack(lev0);

    // Sum the boundaries to recapture any solids moved across
    // grid boundaries during the redistribute
    m_ptr[lev0]->SumBoundary(m_pc->Geom(lev0).periodicity());
    m_ptr[lev0]->FillBoundary(m_pc->Geom(lev0).periodicity());

  } else if (m_redist_type == m_StateRedist) {

    Redistribute_SRD(lev0);

  } else {

  }

}

//
// Redistribute solids volume fraction
//
void MFIXDepositionOp::
Redistribute_MaxPack (int const a_lev)
{
  BL_PROFILE("MFIXDepositionOp::Redistribute_MaxPack");

  const DistributionMapping& dm = m_ptr[a_lev]->DistributionMap();
  const BoxArray&            ba = m_ptr[a_lev]->boxArray();

  int const ncomp = m_ptr[a_lev]->nComp();
  int const ngrow = 1;

  Real const max_eps(m_redist_max_pack);

  MultiFab ptr_copy(ba, dm, ncomp, ngrow);
  MultiFab::Copy(ptr_copy, *m_ptr[a_lev], 0, 0, ncomp, ngrow);

  MultiFab scale_fab(ba, dm, 1, ngrow);
  scale_fab.setVal(0.);

  AMREX_ALWAYS_ASSERT(a_lev == 0);

  const MultiFab* volfrac = &(m_lev0_factory->getVolFrac());
  const FabArray<EBCellFlagFab>* flagFab = &(m_lev0_factory->getMultiEBCellFlagFab());

  for (MFIter mfi(*m_eps[a_lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    // We don't want to do this for ghost cells
    const Box& bx = mfi.tilebox();

    // We are only interested in redistributing excessive particle volume
    // from small cells.

    if ( (*flagFab)[mfi].getType(amrex::grow(bx,1)) != FabType::covered &&
         (*flagFab)[mfi].getType(amrex::grow(bx,1)) != FabType::regular ) {

      // We don't want to loop over ghost cells. Only redistribute
      // solids volume that is locally owned.

      Box domain(m_pc->Geom(a_lev).Domain());
      const amrex::Dim3 dom_low  = amrex::lbound(domain);
      const amrex::Dim3 dom_high = amrex::ubound(domain);

      Array4<const EBCellFlag> const& flags = (*flagFab)[mfi].array();
      Array4<const Real> const& vfrac = volfrac->array(mfi);

      Array4<Real> const& ep_s = m_eps[a_lev]->array(mfi);

      const int cyclic_x = m_pc->Geom(a_lev).isPeriodic(0);
      const int cyclic_y = m_pc->Geom(a_lev).isPeriodic(1);
      const int cyclic_z = m_pc->Geom(a_lev).isPeriodic(2);

      const Box& grow_bx2 = amrex::grow(bx,2);
      IntVect mask_bxg2_lo(grow_bx2.smallEnd());
      IntVect mask_bxg2_hi(grow_bx2.bigEnd());

      if(!cyclic_x) {
        mask_bxg2_lo[0] = amrex::max(mask_bxg2_lo[0], dom_low.x);
        mask_bxg2_hi[0] = amrex::min(mask_bxg2_hi[0], dom_high.x);
      }

      if(!cyclic_y) {
        mask_bxg2_lo[1] = amrex::max(mask_bxg2_lo[1], dom_low.y);
        mask_bxg2_hi[1] = amrex::min(mask_bxg2_hi[1], dom_high.y);
      }

      if(!cyclic_z) {
        mask_bxg2_lo[2] = amrex::max(mask_bxg2_lo[2], dom_low.z);
        mask_bxg2_hi[2] = amrex::min(mask_bxg2_hi[2], dom_high.z);
      }

      // Box "mask_bxg2" is used to restrict were we redistribute the overflow.
      // The following is what we want to do:
      // -- Mask ghost cells when the BCs are not periodic
      // -- Mask cells we are going to redistribute (ep_s > max_eps)
      Box mask_bxg2(mask_bxg2_lo, mask_bxg2_hi);

      Array4<Real> const& mf_redist = m_ptr[a_lev]->array(mfi);

      Array4<Real> const& duplicate = ptr_copy.array(mfi);
      Array4<Real> const& scale_array = scale_fab.array(mfi);

      const Box& grow_bx1 = amrex::grow(bx,1);

      const Box mask_bxg1 = grow_bx1&mask_bxg2;

      amrex::ParallelFor(mask_bxg1,
        [bx, flags,ep_s,mf_redist,vfrac,duplicate,scale_array,max_eps,ncomp,
         mask_bxg2, mask_bxg1] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if(flags(i,j,k).isSingleValued() && ep_s(i,j,k) > max_eps)
        {
          Real sum_vfrac_eps = 0.0;
          Real sum_vfrac     = 0.0;

          for(int ii(-1); ii <= 1; ii++)
          for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++) {
            if((ii != 0 || jj != 0 || kk != 0 ) &&
               flags(i,j,k).isConnected({ii,jj,kk}) &&
               mask_bxg2.contains(IntVect(i+ii,j+jj,k+kk)) &&
               ((!flags(i+ii,j+jj,k+kk).isSingleValued()) ||
                (ep_s(i+ii,j+jj,k+kk) <= max_eps)))
            {
              sum_vfrac     += vfrac(i+ii,j+jj,k+kk);
              sum_vfrac_eps += vfrac(i+ii,j+jj,k+kk)*ep_s(i+ii,j+jj,k+kk);
            }
          }

          // Average volume fraction in the neighborhood around of the packed
          // cell. This value is capped by the user-defined max pack value.
          Real avg_eps = amrex::min(max_eps, sum_vfrac_eps / sum_vfrac);

          // Fraction of material we want to redistribute
          Real scale = amrex::max(0.0, ep_s(i,j,k) - avg_eps) / ep_s(i,j,k);
          scale_array(i,j,k) = scale;

          for(int n(0); n < ncomp; n++) {

            // This is the amount of the multifab we are redistributing
            // that we are moving into the packed cell's neighborhood.
            Real overflow =
              duplicate(i,j,k,n) * scale * vfrac(i,j,k) / sum_vfrac;

            for(int kk(-1); kk <= 1; kk++)
            for(int jj(-1); jj <= 1; jj++)
            for(int ii(-1); ii <= 1; ii++) {
              if((ii != 0 || jj != 0 || kk != 0) &&
                 flags(i,j,k).isConnected({ii,jj,kk}) &&
                 mask_bxg1.contains(IntVect(i+ii,j+jj,k+kk)) &&
                 ((!flags(i+ii,j+jj,k+kk).isSingleValued()) ||
                  (ep_s(i+ii,j+jj,k+kk) <= max_eps)))
              {
                HostDevice::Atomic::Add(&mf_redist(i+ii,j+jj,k+kk,n), overflow);
              }
            }
          }
       }
      });

      Gpu::synchronize();

      amrex::ParallelFor(bx,
        [flags,ep_s,mf_redist,scale_array,max_eps,ncomp]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if(flags(i,j,k).isSingleValued() && ep_s(i,j,k) > max_eps)
        {
          const Real scale = scale_array(i,j,k);

          for(int n(0); n < ncomp; n++) {
            Real redist_val = mf_redist(i,j,k,n);

            // Account for the change in material in the source cell
            Real delta = -scale*redist_val;
            redist_val += delta;
            mf_redist(i,j,k,n) = redist_val;
          }
        }
      });
    }// not regular or covered
  }// MFIter
}

//
// Redistribute solids volume fraction using amrex's SRD algorithm
//
void MFIXDepositionOp::
Redistribute_SRD (int const a_lev)
{
  BL_PROFILE("MFIXDepositionOp::Redistribute_SRD");

  Real const min_vfrac(m_redist_vfrac);

  MultiFab* mf_tmp;
  mf_tmp = (MFHelpers::createFrom(*m_ptr[a_lev])).release();

  const int ncomp = m_ptr[a_lev]->nComp();

  Gpu::DeviceVector<amrex::BCRec> bcs_dummy;
  bcs_dummy.resize(AMREX_SPACEDIM);

  const FabArray<EBCellFlagFab>* flagFab = &(m_lev0_factory->getMultiEBCellFlagFab());

  for (MFIter mfi(*m_ptr[a_lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    Box const& bx = mfi.tilebox();

    Array4<Real> new_arr = mf_tmp->array(mfi);
    Array4<Real> orig_arr = m_ptr[a_lev]->array(mfi);

    bool regular = ((*flagFab)[mfi].getType(amrex::grow(bx,1)) == FabType::regular);
    bool covered = ((*flagFab)[mfi].getType(amrex::grow(bx,1)) == FabType::covered);

    if (!regular && !covered) {

      auto const& flags = (*flagFab)[mfi].const_array();

      auto const& vfrac = (m_lev0_factory->getVolFrac()).const_array(mfi);
      auto const& ccc = (m_lev0_factory->getCentroid()).const_array(mfi);

      auto const& apx = (m_lev0_factory->getAreaFrac()[0])->const_array(mfi);
      auto const& apy = (m_lev0_factory->getAreaFrac()[1])->const_array(mfi);
      auto const& apz = (m_lev0_factory->getAreaFrac()[2])->const_array(mfi);

      auto const& fcx = (m_lev0_factory->getFaceCent()[0])->const_array(mfi);
      auto const& fcy = (m_lev0_factory->getFaceCent()[1])->const_array(mfi);
      auto const& fcz = (m_lev0_factory->getFaceCent()[2])->const_array(mfi);

      ApplyInitialRedistribution(bx, ncomp, new_arr, orig_arr, flags,
        apx, apy, apz, vfrac, fcx, fcy, fcz, ccc, bcs_dummy.data(),
        m_pc->Geom(a_lev), "StateRedist", 2, min_vfrac);

    }
  }

  // Copy back into the original MultiFab
  MultiFab::Copy(*m_ptr[a_lev], *mf_tmp, 0,0, ncomp, 0);

  // we no longer need the copy.
  delete mf_tmp;
}
