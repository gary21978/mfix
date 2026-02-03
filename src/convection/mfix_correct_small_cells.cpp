#include <mfix.H>
#include <mfix_eb_parms.H>

void
mfix::mfix_correct_small_cells (Vector<MultiFab*      > const& vel_in,
                                Vector<MultiFab const*> const& ep_u_mac,
                                Vector<MultiFab const*> const& ep_v_mac,
                                Vector<MultiFab const*> const& ep_w_mac,
                                Vector<MultiFab const*> const& eb_flow_vel)
{
  BL_PROFILE("mfix::mfix_correct_small_cells");

  // Get EB geometric info
  Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
  Real sm_vf = m_correction_small_volfrac;

  for (int lev = 0; lev < nlev(); lev++)
  {
    areafrac = m_eb->factory()[lev]->getAreaFrac();

    for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      const Box bx = mfi.tilebox();

      const EBFArrayBox& vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
      const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

      // Face-centered velocity components
      const auto& umac_fab = (ep_u_mac[lev])->array(mfi);
      const auto& vmac_fab = (ep_v_mac[lev])->array(mfi);
      const auto& wmac_fab = (ep_w_mac[lev])->array(mfi);

      if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
      {
        // do nothing
      }

      // No cut cells in this FAB
      else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
      {
        // do nothing
      }

      // Cut cells in this FAB
      else
      {
        const MultiFab* volfrac = &(m_eb->factory()[lev] -> getVolFrac());

        // Face-centered areas
        const auto& apx_fab = areafrac[0]->array(mfi);
        const auto& apy_fab = areafrac[1]->array(mfi);
        const auto& apz_fab = areafrac[2]->array(mfi);

        const auto& vfrac_fab = volfrac->array(mfi);
        const auto& ccvel_fab = vel_in[lev]->array(mfi);

        if(!m_embedded_boundaries.has_flow()) {

          // This FAB has cut cells -- we define the centroid value in terms
          // of the MAC velocities onfaces
          amrex::ParallelFor(bx,
            [vfrac_fab,apx_fab,apy_fab,apz_fab,ccvel_fab,umac_fab,vmac_fab,wmac_fab,sm_vf]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Real vfrac = vfrac_fab(i,j,k);
            if (vfrac > 0.0 && vfrac < sm_vf)
            {
              const Real apx_mns = apx_fab(i,j,k);
              const Real apx_pls = apx_fab(i+1,j,k);

              const Real apy_mns = apy_fab(i,j,k);
              const Real apy_pls = apy_fab(i,j+1,k);

              const Real apz_mns = apz_fab(i,j,k);
              const Real apz_pls = apz_fab(i,j,k+1);

              ccvel_fab(i,j,k,0) = (apx_mns == 0.0 && apx_pls == 0.0) ? 0.0 :
                                      (apx_mns*umac_fab(i,j,k)+apx_pls*umac_fab(i+1,j,k))/(apx_mns+apx_pls);
              ccvel_fab(i,j,k,1) = (apy_mns == 0.0 && apy_pls == 0.0) ? 0.0 :
                                      (apy_mns*vmac_fab(i,j,k)+apy_pls*vmac_fab(i,j+1,k))/(apy_mns+apy_pls);
              ccvel_fab(i,j,k,2) = (apz_mns == 0.0 && apz_pls == 0.0) ? 0.0:
                                      (apz_mns*wmac_fab(i,j,k)+apz_pls*wmac_fab(i,j,k+1))/(apz_mns+apz_pls);
            }
          });

        } else { // EB has flow

          Array4<Real const> const& eb_vel = (eb_flow_vel[lev])->const_array(mfi);
          Array4<Real const> const& barea = (&(m_eb->factory()[lev]->getBndryArea()))->const_array(mfi);

          // This FAB has cut cells -- we define the centroid value in terms
          // of the MAC velocities onfaces
          amrex::ParallelFor(bx,
            [vfrac_fab,apx_fab,apy_fab,apz_fab,ccvel_fab,umac_fab,vmac_fab,wmac_fab,
             eb_vel,barea,sm_vf]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Real vfrac = vfrac_fab(i,j,k);
            if (vfrac > 0.0 && vfrac < sm_vf)
            {
              const Real apx_mns = apx_fab(i  ,j  ,k  );
              const Real apx_pls = apx_fab(i+1,j  ,k  );
              const Real apy_mns = apy_fab(i  ,j  ,k  );
              const Real apy_pls = apy_fab(i  ,j+1,k  );
              const Real apz_mns = apz_fab(i  ,j  ,k  );
              const Real apz_pls = apz_fab(i  ,j  ,k+1);

              const Real ucc = ccvel_fab(i,j,k,0);
              const Real vcc = ccvel_fab(i,j,k,1);
              const Real wcc = ccvel_fab(i,j,k,2);

              ccvel_fab(i,j,k,0) = (apx_mns == 0.0 && apx_pls == 0.0) ? 0.0 :
                                      (apx_mns*umac_fab(i,j,k)+apx_pls*umac_fab(i+1,j,k))/(apx_mns+apx_pls);
              ccvel_fab(i,j,k,1) = (apy_mns == 0.0 && apy_pls == 0.0) ? 0.0 :
                                      (apy_mns*vmac_fab(i,j,k)+apy_pls*vmac_fab(i,j+1,k))/(apy_mns+apy_pls);
              ccvel_fab(i,j,k,2) = (apz_mns == 0.0 && apz_pls == 0.0) ? 0.0:
                                      (apz_mns*wmac_fab(i,j,k)+apz_pls*wmac_fab(i,j,k+1))/(apz_mns+apz_pls);

              const Real u_eb = eb_vel(i,j,k,0);
              const Real v_eb = eb_vel(i,j,k,1);
              const Real w_eb = eb_vel(i,j,k,2);

              const Real eb_vel_mag = std::sqrt(u_eb*u_eb + v_eb*v_eb + w_eb*w_eb);

              constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

              if (eb_vel_mag > tolerance) {

                 if ( u_eb > tolerance ) {
                     ccvel_fab(i,j,k,0) = amrex::min(ucc, u_eb);
                 } else if( u_eb < -tolerance) {
                     ccvel_fab(i,j,k,0) = amrex::max(ucc, u_eb);
                 }

                 if ( v_eb > tolerance ) {
                     ccvel_fab(i,j,k,1) = amrex::min(vcc, v_eb);
                 } else if( v_eb < -tolerance) {
                     ccvel_fab(i,j,k,1) = amrex::max(vcc, v_eb);
                 }

                 if ( w_eb > tolerance ) {
                     ccvel_fab(i,j,k,2) = amrex::min(wcc, w_eb);
                 } else if( w_eb < -tolerance) {
                     ccvel_fab(i,j,k,2) = amrex::max(wcc, w_eb);
                 }


              }
            }
          });

        }

      }



    }
  }
}
