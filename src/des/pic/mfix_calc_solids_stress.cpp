#include <AMReX_BC_TYPES.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>

#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_pic_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_K.H>
#include <mfix_mf_helpers.H>

void mfix::
MFIX_CalcSolidsStress ( Vector< MultiFab* >& ep_s_in,
                        RealVect& gravity_in,
                        LoadBalance* const a_loadbalance)
{
  BL_PROFILE("mfix::MFIX_CalcSolidsStress()");

  //const Real covered_val = 9.8765e300; UNUSED

  for( int lev(0); lev < nlev(); lev ++)
  {

    const BoxArray&            pba = pc->ParticleBoxArray(lev);
    const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

    // Solids stress gradient. Note that it has an extra component to
    // store a copy of the solids volume fraction for interpolation.
    MultiFab tau;
    tau.define(pba, pdm, 1, 1, MFInfo(), *(m_eb->particle_factory()[lev]));
    tau.setVal(-1.0e100, 0, 1, tau.nGrow());

    const auto& factory = *(m_eb->particle_factory()[lev]);
    const auto& flags = factory.getMultiEBCellFlagFab();

    const Real Ps0 = m_pic.Ps();
    const Real beta = m_pic.beta();
    const Real ep_cp = m_pic.ep_cp();
    const Real small_number = m_pic.small_number();

    constexpr Real max_eps = 0.95;


    for (MFIter mfi(tau, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      const Box& bx = mfi.tilebox();

      Array4<const Real> const& ep_s_arr = ep_s_in[lev]->const_array(mfi);
      Array4<      Real> const& tau_arr   = tau.array(mfi);

      const auto& flagfab = flags[mfi];

      // Fully covered FAB
      if (flagfab.getType(amrex::grow(bx,0)) == FabType::covered ) {
        amrex::ParallelFor(bx,[tau_arr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           tau_arr(i,j,k) = 0.0;
        });

      } else if (flagfab.getType(amrex::grow(bx,1)) == FabType::regular ) {

        amrex::ParallelFor(bx,
          [ep_s_arr,tau_arr,Ps0,beta,ep_cp,small_number,max_eps]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

          const Real bounded_eps = amrex::min(max_eps, ep_s_arr(i,j,k));
          tau_arr(i,j,k) = solids_stress(Ps0, beta, ep_cp, small_number, bounded_eps);

        });

      } else {
        const auto& flagsarr = flagfab.array();

        amrex::ParallelFor(bx,
          [ep_s_arr,tau_arr,Ps0,beta,ep_cp,small_number,max_eps,flagsarr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (flagsarr(i,j,k).isCovered()) {
            tau_arr(i,j,k) = 0.0;

          } else {

            const Real bounded_eps = amrex::min(max_eps, ep_s_arr(i,j,k));
            tau_arr(i,j,k) = solids_stress(Ps0, beta, ep_cp, small_number, bounded_eps);

          }
        });
      }
    }

    // Set all values outside the domain to a crazy value to make
    // catching errors easier to find.
    tau.setDomainBndry(-2.0e100,geom[lev]);

    // Calculate the pressure in flow boundary cells. This is needed
    // for inflow faces that may need to support the weight of the bed.

    Box domain(geom[lev].Domain());

    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    for (MFIter mfi(tau, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      Array4<Real> const& ep_s_arr = ep_s_in[lev]->array(mfi);
      Array4<Real> const& tau_arr = tau.array(mfi);

      IntVect tau_lo(tau[mfi].loVect());
      IntVect tau_hi(tau[mfi].hiVect());

      const Box tau_bx(tau_lo, tau_hi);

      Array4<const int> const& bct_ilo = bc_list.bc_ilo[lev]->array();
      Array4<const int> const& bct_ihi = bc_list.bc_ihi[lev]->array();
      Array4<const int> const& bct_jlo = bc_list.bc_jlo[lev]->array();
      Array4<const int> const& bct_jhi = bc_list.bc_jhi[lev]->array();
      Array4<const int> const& bct_klo = bc_list.bc_klo[lev]->array();
      Array4<const int> const& bct_khi = bc_list.bc_khi[lev]->array();

      const int nlft = amrex::max(0, dom_lo[0]-tau_lo[0]);
      const int nbot = amrex::max(0, dom_lo[1]-tau_lo[1]);
      const int ndwn = amrex::max(0, dom_lo[2]-tau_lo[2]);

      const int nrgt = amrex::max(0, tau_hi[0]-dom_hi[0]);
      const int ntop = amrex::max(0, tau_hi[1]-dom_hi[1]);
      const int nup  = amrex::max(0, tau_hi[2]-dom_hi[2]);

      // Create InVects for following 2D Boxes
      IntVect bx_yz_lo_lo_2D(tau_lo), bx_yz_lo_hi_2D(tau_hi);
      IntVect bx_yz_hi_lo_2D(tau_lo), bx_yz_hi_hi_2D(tau_hi);
      IntVect bx_xz_lo_lo_2D(tau_lo), bx_xz_lo_hi_2D(tau_hi);
      IntVect bx_xz_hi_lo_2D(tau_lo), bx_xz_hi_hi_2D(tau_hi);
      IntVect bx_xy_lo_lo_2D(tau_lo), bx_xy_lo_hi_2D(tau_hi);
      IntVect bx_xy_hi_lo_2D(tau_lo), bx_xy_hi_hi_2D(tau_hi);

      // Fix lo and hi limits
      bx_yz_lo_lo_2D[0] = dom_lo[0]-1;
      bx_yz_lo_hi_2D[0] = dom_lo[0]-1;
      bx_yz_hi_lo_2D[0] = dom_hi[0]+1;
      bx_yz_hi_hi_2D[0] = dom_hi[0]+1;

      bx_xz_lo_lo_2D[1] = dom_lo[1]-1;
      bx_xz_lo_hi_2D[1] = dom_lo[1]-1;
      bx_xz_hi_lo_2D[1] = dom_hi[1]+1;
      bx_xz_hi_hi_2D[1] = dom_hi[1]+1;

      bx_xy_lo_lo_2D[2] = dom_lo[2]-1;
      bx_xy_lo_hi_2D[2] = dom_lo[2]-1;
      bx_xy_hi_lo_2D[2] = dom_hi[2]+1;
      bx_xy_hi_hi_2D[2] = dom_hi[2]+1;

      // Create 2D boxes for GPU loops
      const Box bx_yz_lo_2D(bx_yz_lo_lo_2D, bx_yz_lo_hi_2D);
      const Box bx_yz_hi_2D(bx_yz_hi_lo_2D, bx_yz_hi_hi_2D);

      const Box bx_xz_lo_2D(bx_xz_lo_lo_2D, bx_xz_lo_hi_2D);
      const Box bx_xz_hi_2D(bx_xz_hi_lo_2D, bx_xz_hi_hi_2D);

      const Box bx_xy_lo_2D(bx_xy_lo_lo_2D, bx_xy_lo_hi_2D);
      const Box bx_xy_hi_2D(bx_xy_hi_lo_2D, bx_xy_hi_hi_2D);

      if (nlft > 0) {

        const Real bound_eps = ((gravity_in[0] < 0.0) ? 0.975*ep_cp : 0.0 );

        amrex::ParallelFor(bx_yz_lo_2D,
        [bct_ilo,dom_lo,tau_arr,ep_s_arr,Ps0,beta,ep_cp, small_number,bound_eps]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

          ep_s_arr(i,j,k) = ep_s_arr(dom_lo[0],j,k);

          if((bct == BCList::pinf) || (bct == BCList::minf))
            ep_s_arr(i,j,k) = amrex::max(bound_eps, ep_s_arr(dom_lo[0],j,k));

          tau_arr(i,j,k) = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));

        });
      } // nlft


      if (nrgt > 0) {

        const Real bound_eps = ((gravity_in[0] > 0.0) ? 0.975*ep_cp : 0.0 );

        amrex::ParallelFor(bx_yz_hi_2D,
        [bct_ihi,dom_hi,tau_arr,ep_s_arr, Ps0,beta,ep_cp,small_number,bound_eps]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

            const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

            ep_s_arr(i,j,k) = ep_s_arr(dom_hi[0],j,k);

            if((bct == BCList::pinf) || (bct == BCList::minf))
              ep_s_arr(i,j,k) = amrex::max(bound_eps, ep_s_arr(dom_hi[0],j,k));

            tau_arr(i,j,k) = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));
        });
      } // nrgt


      if (nbot > 0) {

        const Real bound_eps = ((gravity_in[1] < 0.0) ? 0.975*ep_cp : 0.0 );

        amrex::ParallelFor(bx_xz_lo_2D,
        [bct_jlo,dom_lo,tau_arr,ep_s_arr,Ps0,beta,ep_cp,small_number,bound_eps]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

          ep_s_arr(i,j,k) = ep_s_arr(i,dom_lo[1],k);

          if((bct == BCList::pinf) || (bct == BCList::minf))
            ep_s_arr(i,j,k) = amrex::max(bound_eps, ep_s_arr(i,dom_lo[1],k));

          tau_arr(i,j,k) = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));

        });
      } // nbot


      if (ntop > 0) {

        const Real bound_eps = ((gravity_in[1] > 0.0) ? 0.975*ep_cp : 0.0 );

        amrex::ParallelFor(bx_xz_hi_2D,
        [bct_jhi,dom_hi,tau_arr,ep_s_arr, Ps0,beta,ep_cp,small_number,bound_eps]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

          ep_s_arr(i,j,k) = ep_s_arr(i,dom_hi[1],k);

          if((bct == BCList::pinf) || (bct == BCList::minf))
            ep_s_arr(i,j,k) = amrex::max(bound_eps, ep_s_arr(i,dom_hi[1],k));

          tau_arr(i,j,k) = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));
        });
      } // ntop


      if (ndwn > 0) {

        const Real bound_eps = ((gravity_in[2] < 0.0) ? 0.975*ep_cp : 0.0 );

        amrex::ParallelFor(bx_xy_lo_2D,
        [bct_klo,dom_lo,tau_arr,ep_s_arr, Ps0,beta,ep_cp, small_number,bound_eps]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_klo(i,j,dom_lo[2]-1,0);

          ep_s_arr(i,j,k) = ep_s_arr(i,j,dom_lo[2]);

          if((bct == BCList::pinf) || (bct == BCList::minf))
            ep_s_arr(i,j,k) = amrex::max(bound_eps, ep_s_arr(i,j,dom_lo[2]));

          tau_arr(i,j,k) = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));

        });
      } // ndwn


      if (nup  > 0) {

        const Real bound_eps = ((gravity_in[2] > 0.0) ? 0.975*ep_cp : 0.0 );

        amrex::ParallelFor(bx_xy_hi_2D,
        [bct_khi,dom_hi,tau_arr,ep_s_arr, Ps0,beta,ep_cp,small_number,bound_eps]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_khi(i,j,dom_hi[2]+1,0);

          ep_s_arr(i,j,k) = ep_s_arr(i,j,dom_hi[2]);
          if((bct == BCList::pinf) || (bct == BCList::minf))
            ep_s_arr(i,j,k) = amrex::max(bound_eps, ep_s_arr(i,j,dom_hi[2]));

          tau_arr(i,j,k) = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));
        });
      } // nup

    }

    tau.FillBoundary(geom[lev].periodicity());


// Some debugging code
#if 0
    for (MFIter mfi(tau, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      const Box& bx = mfi.validbox();

      if(bx.contains(epg_cell)) {

        for(int ii(-1); ii<=1; ii++){
            for(int kk(-1); kk<=1; kk++){
              for(int jj(-1); jj<=1; jj++){

              const int i = epg_cell[0] + ii;
              const int j = epg_cell[1] + jj;
              const int k = epg_cell[2] + kk;

              IntVect ijk = {i,j,k};
              //value = mf[mfi](iv,icomp); // icomp : component index
              if(amrex::grow(bx,2).contains(ijk)) {
                        amrex::Print(Print::AllProcs) << ijk << ", "
                          << ep_s_in[lev]->const_array(mfi)(ijk,0) << ", "
                          << tau.const_array(mfi)(ijk,0) << ", "
                          << "\n";
                }

              }
            }
          }
      }
    }
#endif


    // We now have the solids stress on the field. We now need to interpolate
    // this to the parcel positions. We are going to store this for the
    // parcel in "omega" in the particle container.

    using MFIXParIter = MFIXParticleContainer::MFIXParIter;

    {

      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const RealVect plo(plo_array[0], plo_array[1], plo_array[2]);


      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        // Timer used for load-balancing
        Real const timer_start( ParallelDescriptor::second() );

        auto& particles = pti.GetArrayOfStructs();
        MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();

        auto& ptile = pti.GetParticleTile();
        auto ptile_data = ptile.getParticleTileData();
        auto runtime_idxs = pc->m_runtimeRealData;

        const int np = particles.size();

        Box bx = pti.tilebox ();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  ep_s_fab = static_cast<EBFArrayBox const&>((*ep_s_in[lev])[pti]);
        const EBCellFlagFab&  flags_loc = ep_s_fab.getEBCellFlagFab();

        Array4<const Real> const& ep_s_arr = ep_s_in[lev]->const_array(pti);
        Array4<const Real> const& tau_arr   = tau.array(pti);
        const auto& flags_array = flags_loc.array();

        if (flags_loc.getType(amrex::grow(bx,1)) == FabType::covered)
        {

          // We shouldn't have this case but if we do -- zero the stress.
          amrex::ParallelFor(np, [pstruct, p_realarray, ptile_data, runtime_idxs]
            AMREX_GPU_DEVICE (int pid) noexcept
            {
              ptile_data.m_runtime_rdata[runtime_idxs.ep_s][pid] = 0.0;

              p_realarray[SoArealData::omegax][pid] = 0.0;
              p_realarray[SoArealData::omegay][pid] = 0.0;
              p_realarray[SoArealData::omegaz][pid] = 0.0;
            });
        }
        else if (flags_loc.getType(amrex::grow(bx,1)) == FabType::regular)
        {
          amrex::ParallelFor(np, [pstruct, p_realarray, ptile_data, runtime_idxs,
              ep_s_arr, tau_arr, plo, dxi]
            AMREX_GPU_DEVICE (int pid) noexcept
          {

            MFIXParticleContainer::ParticleType& particle = pstruct[pid];

            const Real lx = (particle.pos(0) - plo[0]) * dxi[0] + 0.5;
            const Real ly = (particle.pos(1) - plo[1]) * dxi[1] + 0.5;
            const Real lz = (particle.pos(2) - plo[2]) * dxi[2] + 0.5;

            const int i = static_cast<int>(amrex::Math::floor(lx));
            const int j = static_cast<int>(amrex::Math::floor(ly));
            const int k = static_cast<int>(amrex::Math::floor(lz));

            const Real wx_hi(lx - static_cast<Real>(i));
            const Real wy_hi(ly - static_cast<Real>(j));
            const Real wz_hi(lz - static_cast<Real>(k));

            const Real wx_lo(1.0 - wx_hi);
            const Real wy_lo(1.0 - wy_hi);
            const Real wz_lo(1.0 - wz_hi);

            AMREX_ASSERT(ep_s_arr(i-1,j-1, k-1) >= 0.);
            AMREX_ASSERT(ep_s_arr(i-1,j-1, k  ) >= 0.);
            AMREX_ASSERT(ep_s_arr(i-1,j  , k-1) >= 0.);
            AMREX_ASSERT(ep_s_arr(i-1,j  , k  ) >= 0.);
            AMREX_ASSERT(ep_s_arr(i  ,j-1, k-1) >= 0.);
            AMREX_ASSERT(ep_s_arr(i  ,j-1, k  ) >= 0.);
            AMREX_ASSERT(ep_s_arr(i  ,j  , k-1) >= 0.);
            AMREX_ASSERT(ep_s_arr(i  ,j  , k  ) >= 0.);

            AMREX_ASSERT(tau_arr(i-1,j-1, k-1) >= 0.);
            AMREX_ASSERT(tau_arr(i-1,j-1, k  ) >= 0.);
            AMREX_ASSERT(tau_arr(i-1,j  , k-1) >= 0.);
            AMREX_ASSERT(tau_arr(i-1,j  , k  ) >= 0.);
            AMREX_ASSERT(tau_arr(i  ,j-1, k-1) >= 0.);
            AMREX_ASSERT(tau_arr(i  ,j-1, k  ) >= 0.);
            AMREX_ASSERT(tau_arr(i  ,j  , k-1) >= 0.);
            AMREX_ASSERT(tau_arr(i  ,j  , k  ) >= 0.);

            const Real ep_s_loc =
              + wx_lo*wy_lo*wz_lo*ep_s_arr(i-1,j-1, k-1)
              + wx_lo*wy_lo*wz_hi*ep_s_arr(i-1,j-1, k  )
              + wx_lo*wy_hi*wz_lo*ep_s_arr(i-1,j  , k-1)
              + wx_lo*wy_hi*wz_hi*ep_s_arr(i-1,j  , k  )
              + wx_hi*wy_lo*wz_lo*ep_s_arr(i  ,j-1, k-1)
              + wx_hi*wy_lo*wz_hi*ep_s_arr(i  ,j-1, k  )
              + wx_hi*wy_hi*wz_lo*ep_s_arr(i  ,j  , k-1)
              + wx_hi*wy_hi*wz_hi*ep_s_arr(i  ,j  , k  );

            const Real dtaudx = dxi[0]*(
              - wy_lo*wz_lo*tau_arr(i-1,j-1, k-1)
              - wy_lo*wz_hi*tau_arr(i-1,j-1, k  )
              - wy_hi*wz_lo*tau_arr(i-1,j  , k-1)
              - wy_hi*wz_hi*tau_arr(i-1,j  , k  )
              + wy_lo*wz_lo*tau_arr(i  ,j-1, k-1)
              + wy_lo*wz_hi*tau_arr(i  ,j-1, k  )
              + wy_hi*wz_lo*tau_arr(i  ,j  , k-1)
              + wy_hi*wz_hi*tau_arr(i  ,j  , k  ));

            const Real dtaudy = dxi[1]*(
              - wx_lo*wz_lo*tau_arr(i-1,j-1, k-1)
              - wx_lo*wz_hi*tau_arr(i-1,j-1, k  )
              + wx_lo*wz_lo*tau_arr(i-1,j  , k-1)
              + wx_lo*wz_hi*tau_arr(i-1,j  , k  )
              - wx_hi*wz_lo*tau_arr(i  ,j-1, k-1)
              - wx_hi*wz_hi*tau_arr(i  ,j-1, k  )
              + wx_hi*wz_lo*tau_arr(i  ,j  , k-1)
              + wx_hi*wz_hi*tau_arr(i  ,j  , k  ));

            const Real dtaudz = dxi[2]*(
              - wx_lo*wy_lo*tau_arr(i-1,j-1, k-1)
              + wx_lo*wy_lo*tau_arr(i-1,j-1, k  )
              - wx_lo*wy_hi*tau_arr(i-1,j  , k-1)
              + wx_lo*wy_hi*tau_arr(i-1,j  , k  )
              - wx_hi*wy_lo*tau_arr(i  ,j-1, k-1)
              + wx_hi*wy_lo*tau_arr(i  ,j-1, k  )
              - wx_hi*wy_hi*tau_arr(i  ,j  , k-1)
              + wx_hi*wy_hi*tau_arr(i  ,j  , k  ));


            // Store local solids volume fraction in unused array location.
            ptile_data.m_runtime_rdata[runtime_idxs.ep_s][pid] = ep_s_loc;

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_loc >= 0.,"Check valid ep_s_loc");

            p_realarray[SoArealData::omegax][pid] = dtaudx;
            p_realarray[SoArealData::omegay][pid] = dtaudy;
            p_realarray[SoArealData::omegaz][pid] = dtaudz;
          });
        }
        else // FAB not all regular
        {
          amrex::ParallelFor(np,
           [pstruct, p_realarray, ptile_data, runtime_idxs, flags_array, ep_s_arr,
            tau_arr, plo, dxi, Ps0, beta, ep_cp, small_number]
          AMREX_GPU_DEVICE (int pid) noexcept
          {

            MFIXParticleContainer::ParticleType& particle = pstruct[pid];

            const int ip = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
            const int jp = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
            const int kp = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

            if(flags_array(ip,jp,kp).isCovered())
            {
              p_realarray[SoArealData::omegax][pid] = 0.0;
              p_realarray[SoArealData::omegay][pid] = 0.0;
              p_realarray[SoArealData::omegaz][pid] = 0.0;
            // Cut or regular cell and none of the cells in the stencil is covered
            // (Note we can't assume regular cell has no covered cells in the stencil
            //      because of the diagonal case)
            } else {

              const Real lx = (particle.pos(0) - plo[0]) * dxi[0] + 0.5;
              const Real ly = (particle.pos(1) - plo[1]) * dxi[1] + 0.5;
              const Real lz = (particle.pos(2) - plo[2]) * dxi[2] + 0.5;

              // Upper cell in trilinear stencil
              const int i = static_cast<int>(amrex::Math::floor(lx));
              const int j = static_cast<int>(amrex::Math::floor(ly));
              const int k = static_cast<int>(amrex::Math::floor(lz));

              const Real wx_hi(lx - static_cast<Real>(i));
              const Real wy_hi(ly - static_cast<Real>(j));
              const Real wz_hi(lz - static_cast<Real>(k));

              const Real wx_lo(1.0 - wx_hi);
              const Real wy_lo(1.0 - wy_hi);
              const Real wz_lo(1.0 - wz_hi);

              const int di = i - ip; // 0 or 1
              const int dj = j - jp; // 0 or 1
              const int dk = k - kp; // 0 or 1

              Real ep_s_loc;

              Real dtaudx;
              Real dtaudy;
              Real dtaudz;

              Real ep_s_wt = 0.0;

              // Check that all the cells in the stencil are connected to the
              // that contains the parcels centroid. We don't need to check
              // cell (ip,jp,kp) as we've already tested if it is covered.

              if (flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk-1) &&
                  flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk-1) &&
                  flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk-1) &&
                  flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk-1) &&
                  flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk  ) &&
                  flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk  ) &&
                  flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk  ) &&
                  flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk  )) {

                // After sufficient testing, these should be changed to
                // standard AMREX_ASSERT checks
                AMREX_ASSERT(ep_s_arr(i-1,j-1, k-1) >= 0. && ep_s_arr(i-1,j-1, k-1) <  1.);
                AMREX_ASSERT(ep_s_arr(i-1,j-1, k  ) >= 0. && ep_s_arr(i-1,j-1, k  ) <  1.);
                AMREX_ASSERT(ep_s_arr(i-1,j  , k-1) >= 0. && ep_s_arr(i-1,j  , k-1) <  1.);
                AMREX_ASSERT(ep_s_arr(i-1,j  , k  ) >= 0. && ep_s_arr(i-1,j  , k  ) <  1.);
                AMREX_ASSERT(ep_s_arr(i  ,j-1, k-1) >= 0. && ep_s_arr(i  ,j-1, k-1) <  1.);
                AMREX_ASSERT(ep_s_arr(i  ,j-1, k  ) >= 0. && ep_s_arr(i  ,j-1, k  ) <  1.);
                AMREX_ASSERT(ep_s_arr(i  ,j  , k-1) >= 0. && ep_s_arr(i  ,j  , k-1) <  1.);
                AMREX_ASSERT(ep_s_arr(i  ,j  , k  ) >= 0. && ep_s_arr(i  ,j  , k  ) <  1.);

                AMREX_ASSERT(tau_arr(i-1,j-1, k-1) >= 0.);
                AMREX_ASSERT(tau_arr(i-1,j-1, k  ) >= 0.);
                AMREX_ASSERT(tau_arr(i-1,j  , k-1) >= 0.);
                AMREX_ASSERT(tau_arr(i-1,j  , k  ) >= 0.);
                AMREX_ASSERT(tau_arr(i  ,j-1, k-1) >= 0.);
                AMREX_ASSERT(tau_arr(i  ,j-1, k  ) >= 0.);
                AMREX_ASSERT(tau_arr(i  ,j  , k-1) >= 0.);
                AMREX_ASSERT(tau_arr(i  ,j  , k  ) >= 0.);

                ep_s_loc =
                  + wx_lo*wy_lo*wz_lo*ep_s_arr(i-1,j-1, k-1)
                  + wx_lo*wy_lo*wz_hi*ep_s_arr(i-1,j-1, k  )
                  + wx_lo*wy_hi*wz_lo*ep_s_arr(i-1,j  , k-1)
                  + wx_lo*wy_hi*wz_hi*ep_s_arr(i-1,j  , k  )
                  + wx_hi*wy_lo*wz_lo*ep_s_arr(i  ,j-1, k-1)
                  + wx_hi*wy_lo*wz_hi*ep_s_arr(i  ,j-1, k  )
                  + wx_hi*wy_hi*wz_lo*ep_s_arr(i  ,j  , k-1)
                  + wx_hi*wy_hi*wz_hi*ep_s_arr(i  ,j  , k  );

                dtaudx = dxi[0]*(
                  - wy_lo*wz_lo*tau_arr(i-1,j-1, k-1)
                  - wy_lo*wz_hi*tau_arr(i-1,j-1, k  )
                  - wy_hi*wz_lo*tau_arr(i-1,j  , k-1)
                  - wy_hi*wz_hi*tau_arr(i-1,j  , k  )
                  + wy_lo*wz_lo*tau_arr(i  ,j-1, k-1)
                  + wy_lo*wz_hi*tau_arr(i  ,j-1, k  )
                  + wy_hi*wz_lo*tau_arr(i  ,j  , k-1)
                  + wy_hi*wz_hi*tau_arr(i  ,j  , k  ));

                dtaudy = dxi[1]*(
                  - wx_lo*wz_lo*tau_arr(i-1,j-1, k-1)
                  - wx_lo*wz_hi*tau_arr(i-1,j-1, k  )
                  + wx_lo*wz_lo*tau_arr(i-1,j  , k-1)
                  + wx_lo*wz_hi*tau_arr(i-1,j  , k  )
                  - wx_hi*wz_lo*tau_arr(i  ,j-1, k-1)
                  - wx_hi*wz_hi*tau_arr(i  ,j-1, k  )
                  + wx_hi*wz_lo*tau_arr(i  ,j  , k-1)
                  + wx_hi*wz_hi*tau_arr(i  ,j  , k  ));

                dtaudz = dxi[2]*(
                  - wx_lo*wy_lo*tau_arr(i-1,j-1, k-1)
                  + wx_lo*wy_lo*tau_arr(i-1,j-1, k  )
                  - wx_lo*wy_hi*tau_arr(i-1,j  , k-1)
                  + wx_lo*wy_hi*tau_arr(i-1,j  , k  )
                  - wx_hi*wy_lo*tau_arr(i  ,j-1, k-1)
                  + wx_hi*wy_lo*tau_arr(i  ,j-1, k  )
                  - wx_hi*wy_hi*tau_arr(i  ,j  , k-1)
                  + wx_hi*wy_hi*tau_arr(i  ,j  , k  ));

              // At least one of the cells in the stencil is not connected
              } else {

                ep_s_loc = 0.0;

                { ///////////////////////////////// Calculate dtau_dx

                  Real eps_arr[2][2];
                  Real dtau_arr[2][2];
                  int lci[2][2];

                  //const int ii = i + di - 1;

                  for(int kk=0; kk<=1; kk++){
                    for(int jj=0; jj<=1; jj++){

                      eps_arr[1-jj][1-kk] = 0.0;
                      dtau_arr[1-jj][1-kk] = 0.0;
                      lci[1-jj][1-kk]=0;

                      // Both cells have valid values so calculate dtau as usual.
                      if ((!flags_array(i  ,j-jj,k-kk).isCovered()) &&
                          (!flags_array(i-1,j-jj,k-kk).isCovered())){

                        eps_arr[1-jj][1-kk] = wx_hi*ep_s_arr(i  ,j-jj,k-kk) + wx_lo*ep_s_arr(i-1,j-jj,k-kk);
                        dtau_arr[1-jj][1-kk] = (tau_arr(i  ,j-jj,k-kk) - tau_arr(i-1,j-jj,k-kk));
                        lci[1-jj][1-kk]=1;

                      // Particle has index i -- extrap in x to i-1 if needed
                      } else if (di == 0 && !flags_array(i  ,j-jj,k-kk).isCovered()){

                        AMREX_ASSERT( !flags_array(i+1,j-jj,k-kk).isCovered());

                        const Real ep_s_extrap = ep_s_arr(i  ,j-jj,k-kk);

                        const Real tau_lo = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                        eps_arr[1-jj][1-kk] = wx_hi*ep_s_arr(i  ,j-jj,k-kk) + wx_lo*ep_s_extrap;
                        dtau_arr[1-jj][1-kk] = (tau_arr(i  ,j-jj,k-kk) - tau_lo);
                        lci[1-jj][1-kk]=1;

                      // particle has index i-1 -- extrap in x to i if needed
                      } else if (di == 1 && !flags_array(i-1,j-jj,k-kk).isCovered()){

                        AMREX_ASSERT( !flags_array(i-2,j-jj,k-kk).isCovered());

                        const Real ep_s_extrap = ep_s_arr(i-1,j-jj,k-kk);

                        const Real tau_hi = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                        eps_arr[1-jj][1-kk] = wx_hi*ep_s_extrap + wx_lo*ep_s_arr(i-1,j-jj,k-kk);
                        dtau_arr[1-jj][1-kk] = (tau_hi - tau_arr(i-1,j-jj,k-kk));
                        lci[1-jj][1-kk]=1;

                      }

                    } // loop over jj
                  } // loop over kk

                  // How many valid calculations of dtau_dx do we have?
                  const int count = lci[0][0] + lci[1][0] + lci[0][1] + lci[1][1];

                  const Real mij_00 = static_cast<Real>(lci[0][0]);
                  const Real mij_01 = static_cast<Real>(lci[0][1]);
                  const Real mij_10 = static_cast<Real>(lci[1][0]);
                  const Real mij_11 = static_cast<Real>(lci[1][1]);

                  ep_s_wt += static_cast<Real>(count);

                  if ( count == 4 ) {

                    ep_s_loc += 4.0*(
                      + wy_lo*wz_lo*eps_arr[0][0]
                      + wy_lo*wz_hi*eps_arr[0][1]
                      + wy_hi*wz_lo*eps_arr[1][0]
                      + wy_hi*wz_hi*eps_arr[1][1]);

                    dtaudx = dxi[0]*(
                      + wy_lo*wz_lo*dtau_arr[0][0]
                      + wy_lo*wz_hi*dtau_arr[0][1]
                      + wy_hi*wz_lo*dtau_arr[1][0]
                      + wy_hi*wz_hi*dtau_arr[1][1]);

                  } else if ( count == 3 ) {

                    const Real djk_00 = wy_hi*wy_hi + wz_hi*wz_hi;
                    const Real djk_01 = wy_hi*wy_hi + wz_lo*wz_lo;
                    const Real djk_10 = wy_lo*wy_lo + wz_hi*wz_hi;
                    const Real djk_11 = wy_lo*wy_lo + wz_lo*wz_lo;

                    const Real wjk_00 = mij_00 * djk_01 * djk_10 * djk_11;
                    const Real wjk_01 = djk_00 * mij_01 * djk_10 * djk_11;
                    const Real wjk_10 = djk_00 * djk_01 * mij_10 * djk_11;
                    const Real wjk_11 = djk_00 * djk_01 * djk_10 * mij_11;

                    const Real inv_denom = 1.0 / (wjk_00 + wjk_01 + wjk_10 + wjk_11);

                    ep_s_loc += 3.0*(
                      + wjk_00*eps_arr[0][0]
                      + wjk_01*eps_arr[0][1]
                      + wjk_10*eps_arr[1][0]
                      + wjk_11*eps_arr[1][1])*inv_denom;

                    dtaudx = dxi[0]*(
                      + wjk_00*dtau_arr[0][0]
                      + wjk_01*dtau_arr[0][1]
                      + wjk_10*dtau_arr[1][0]
                      + wjk_11*dtau_arr[1][1])*inv_denom;

                  } else if ( count == 2 ) {

                    ep_s_loc += 2.0 * (
                      + mij_00*mij_10*(wy_lo*eps_arr[0][0] + wy_hi*eps_arr[1][0])
                      + mij_00*mij_01*(wz_lo*eps_arr[0][0] + wz_hi*eps_arr[0][1])
                      + mij_11*mij_10*(wz_lo*eps_arr[1][0] + wz_hi*eps_arr[1][1])
                      + mij_11*mij_01*(wy_lo*eps_arr[0][1] + wy_hi*eps_arr[1][1]));

                    dtaudx = dxi[0]*(
                      + mij_00*mij_10*(wy_lo*dtau_arr[0][0] + wy_hi*dtau_arr[1][0])
                      + mij_00*mij_01*(wz_lo*dtau_arr[0][0] + wz_hi*dtau_arr[0][1])
                      + mij_11*mij_10*(wz_lo*dtau_arr[1][0] + wz_hi*dtau_arr[1][1])
                      + mij_11*mij_01*(wy_lo*dtau_arr[0][1] + wy_hi*dtau_arr[1][1]));


                  } else {

                    AMREX_ALWAYS_ASSERT( count == 1);

                    ep_s_loc +=
                      + mij_00*eps_arr[0][0]
                      + mij_01*eps_arr[0][1]
                      + mij_10*eps_arr[1][0]
                      + mij_11*eps_arr[1][1];

                    dtaudx = dxi[0]*(
                      + mij_00*dtau_arr[0][0]
                      + mij_01*dtau_arr[0][1]
                      + mij_10*dtau_arr[1][0]
                      + mij_11*dtau_arr[1][1]);
                  }

                } // END dtau_dx


                { ///////////////////////////////// Calculate dtau_dy

                  Real eps_arr[2][2];
                  Real dtau_arr[2][2];
                  int lci[2][2];

                  //onst int jj = j + dj - 1;

                  for(int kk=0; kk<=1; kk++){
                    for(int ii=0; ii<=1; ii++){

                      eps_arr[1-ii][1-kk] = 0.0;
                      dtau_arr[1-ii][1-kk] = 0.0;
                      lci[1-ii][1-kk]=0;

                      // Both cells have valid values so calculate dtau as usual.
                      if ((!flags_array(i-ii,j  ,k-kk).isCovered()) &&
                          (!flags_array(i-ii,j-1,k-kk).isCovered())){

                        eps_arr[1-ii][1-kk] = wy_hi*ep_s_arr(i-ii,j  ,k-kk) + wy_lo*ep_s_arr(i-ii,j-1,k-kk);
                        dtau_arr[1-ii][1-kk] = (tau_arr(i-ii,j  ,k-kk) - tau_arr(i-ii,j-1,k-kk));
                        lci[1-ii][1-kk]=1;

                      // Particle has index j -- extrap in y to j-1 if needed
                      } else if (dj == 0 && !flags_array(i-ii,j  ,k-kk).isCovered()){

                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( !flags_array(i-ii,j+1,k-kk).isCovered(),"Check y-extrap A");

                        const Real ep_s_extrap = ep_s_arr(i-ii,j  ,k-kk);

                        const Real tau_lo = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                        eps_arr[1-ii][1-kk] = wy_hi*ep_s_arr(i-ii,j  ,k-kk) + wy_lo*ep_s_extrap;
                        dtau_arr[1-ii][1-kk] = (tau_arr(i-ii,j  ,k-kk) - tau_lo);
                        lci[1-ii][1-kk]=1;

                      // particle has index j-1 -- extrap in y to j if needed
                      } else if (dj == 1 && !flags_array(i-ii,j-1,k-kk).isCovered()){

                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( !flags_array(i-ii,j-2,k-kk).isCovered(),"Check y-extrap B");

                        const Real ep_s_extrap = ep_s_arr(i-ii,j-1,k-kk);

                        const Real tau_hi = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                        eps_arr[1-ii][1-kk] = wy_hi*ep_s_extrap + wy_lo*ep_s_arr(i-ii,j-1,k-kk);
                        dtau_arr[1-ii][1-kk] = (tau_hi - tau_arr(i-ii,j-1,k-kk));
                        lci[1-ii][1-kk]=1;

                      }

                    } // loop over ii
                  } // loop over kk

                  // How many valid calculations of dtau_dy do we have?
                  const int count = lci[0][0] + lci[1][0] + lci[0][1] + lci[1][1];

                  const Real mij_00 = static_cast<Real>(lci[0][0]);
                  const Real mij_01 = static_cast<Real>(lci[0][1]);
                  const Real mij_10 = static_cast<Real>(lci[1][0]);
                  const Real mij_11 = static_cast<Real>(lci[1][1]);

                  ep_s_wt += static_cast<Real>(count);

                  if ( count == 4 ) {

                    ep_s_loc += 4.0*(
                      + wx_lo*wz_lo*eps_arr[0][0]
                      + wx_lo*wz_hi*eps_arr[0][1]
                      + wx_hi*wz_lo*eps_arr[1][0]
                      + wx_hi*wz_hi*eps_arr[1][1]);

                    dtaudy = dxi[1]*(
                      + wx_lo*wz_lo*dtau_arr[0][0]
                      + wx_lo*wz_hi*dtau_arr[0][1]
                      + wx_hi*wz_lo*dtau_arr[1][0]
                      + wx_hi*wz_hi*dtau_arr[1][1]);

                  } else if ( count == 3 ) {

                    const Real djk_00 = wx_hi*wx_hi + wz_hi*wz_hi;
                    const Real djk_01 = wx_hi*wx_hi + wz_lo*wz_lo;
                    const Real djk_10 = wx_lo*wx_lo + wz_hi*wz_hi;
                    const Real djk_11 = wx_lo*wx_lo + wz_lo*wz_lo;

                    const Real wjk_00 = mij_00 * djk_01 * djk_10 * djk_11;
                    const Real wjk_01 = djk_00 * mij_01 * djk_10 * djk_11;
                    const Real wjk_10 = djk_00 * djk_01 * mij_10 * djk_11;
                    const Real wjk_11 = djk_00 * djk_01 * djk_10 * mij_11;

                    const Real inv_denom = 1.0 / (wjk_00 + wjk_01 + wjk_10 + wjk_11);

                    ep_s_loc += 3.0*(
                      + wjk_00*eps_arr[0][0]
                      + wjk_01*eps_arr[0][1]
                      + wjk_10*eps_arr[1][0]
                      + wjk_11*eps_arr[1][1])*inv_denom;

                    dtaudy = dxi[1]*(
                      + wjk_00*dtau_arr[0][0]
                      + wjk_01*dtau_arr[0][1]
                      + wjk_10*dtau_arr[1][0]
                      + wjk_11*dtau_arr[1][1])*inv_denom;

                  } else if ( count == 2 ) {

                    ep_s_loc += 2.0 * (
                      + mij_00*mij_10*(wx_lo*eps_arr[0][0] + wx_hi*eps_arr[1][0])
                      + mij_00*mij_01*(wz_lo*eps_arr[0][0] + wz_hi*eps_arr[0][1])
                      + mij_11*mij_10*(wz_lo*eps_arr[1][0] + wz_hi*eps_arr[1][1])
                      + mij_11*mij_01*(wx_lo*eps_arr[0][1] + wx_hi*eps_arr[1][1]));

                    dtaudy = dxi[1]*(
                      + mij_00*mij_10*(wx_lo*dtau_arr[0][0] + wx_hi*dtau_arr[1][0])
                      + mij_00*mij_01*(wz_lo*dtau_arr[0][0] + wz_hi*dtau_arr[0][1])
                      + mij_11*mij_10*(wz_lo*dtau_arr[1][0] + wz_hi*dtau_arr[1][1])
                      + mij_11*mij_01*(wx_lo*dtau_arr[0][1] + wx_hi*dtau_arr[1][1]));

                  } else {

                    AMREX_ALWAYS_ASSERT( count == 1);

                    ep_s_loc +=
                      + mij_00*eps_arr[0][0]
                      + mij_01*eps_arr[0][1]
                      + mij_10*eps_arr[1][0]
                      + mij_11*eps_arr[1][1];

                    dtaudy = dxi[1]*(
                      + mij_00*dtau_arr[0][0]
                      + mij_01*dtau_arr[0][1]
                      + mij_10*dtau_arr[1][0]
                      + mij_11*dtau_arr[1][1]);
                  }

                } // END dtau_dy

                { ///////////////////////////////// Calculate dtau_dz

                  Real eps_arr[2][2];
                  Real dtau_arr[2][2];
                  int lci[2][2];

                  //const int kk = k + dk - 1;

                  for(int jj=0; jj<=1; jj++){
                    for(int ii=0; ii<=1; ii++){

                      eps_arr[1-ii][1-jj] = 0.0;
                      dtau_arr[1-ii][1-jj] = 0.0;
                      lci[1-ii][1-jj]=0;

                      // Both cells have valid values so calculate dtau as usual.
                      if ((!flags_array(i-ii,j-jj,k  ).isCovered()) &&
                          (!flags_array(i-ii,j-jj,k-1).isCovered())){

                        eps_arr[1-ii][1-jj] = wz_hi*ep_s_arr(i-ii,j-jj,k  ) + wz_lo*ep_s_arr(i-ii,j-jj,k-1);
                        dtau_arr[1-ii][1-jj] = (tau_arr(i-ii,j-jj,k  ) - tau_arr(i-ii,j-jj,k-1));
                        lci[1-ii][1-jj]=1;

                      // Particle has index k -- extrap in z to k-1 if needed
                      } else if (dk == 0 && !flags_array(i-ii,j-jj,k  ).isCovered()){

                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( !flags_array(i-ii,j-jj,k+1).isCovered(),"Check z-extrap A");

                        const Real ep_s_extrap = ep_s_arr(i-ii,j-jj,k  );

                        const Real tau_lo = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                        eps_arr[1-ii][1-jj] = wz_hi*ep_s_arr(i-ii,j-jj,k  ) + wz_lo*ep_s_extrap;
                        dtau_arr[1-ii][1-jj] = (tau_arr(i-ii,j-jj,k  ) - tau_lo);
                        lci[1-ii][1-jj]=1;

                      // particle has index k-1 -- extrap in z to k if needed
                      } else if (dk == 1 && !flags_array(i-ii,j-jj,k-1).isCovered()){

                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( !flags_array(i-ii,j-jj,k-2).isCovered(),"Check z-extrap B");

                        const Real ep_s_extrap = ep_s_arr(i-ii,j-jj,k-1);

                        const Real tau_hi = solids_stress(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                        eps_arr[1-ii][1-jj] = wz_hi*ep_s_extrap + wz_lo*ep_s_arr(i-ii,j-jj,k-1);
                        dtau_arr[1-ii][1-jj] = (tau_hi - tau_arr(i-ii,j-jj,k-1));
                        lci[1-ii][1-jj]=1;

                      }

                    } // loop over ii
                  } // loop over jj

                  // How many valid calculations of dtau_dz do we have?
                  const int count = lci[0][0] + lci[1][0] + lci[0][1] + lci[1][1];

                  const Real mij_00 = static_cast<Real>(lci[0][0]);
                  const Real mij_01 = static_cast<Real>(lci[0][1]);
                  const Real mij_10 = static_cast<Real>(lci[1][0]);
                  const Real mij_11 = static_cast<Real>(lci[1][1]);

                  ep_s_wt += static_cast<Real>(count);

                  if ( count == 4 ) {

                    ep_s_loc += 4.0*(
                      + wx_lo*wy_lo*eps_arr[0][0]
                      + wx_lo*wy_hi*eps_arr[0][1]
                      + wx_hi*wy_lo*eps_arr[1][0]
                      + wx_hi*wy_hi*eps_arr[1][1]);

                    dtaudz = dxi[2]*(
                      + wx_lo*wy_lo*dtau_arr[0][0]
                      + wx_lo*wy_hi*dtau_arr[0][1]
                      + wx_hi*wy_lo*dtau_arr[1][0]
                      + wx_hi*wy_hi*dtau_arr[1][1]);

                  } else if ( count == 3 ) {

                    const Real djk_00 = wx_hi*wx_hi + wy_hi*wy_hi;
                    const Real djk_01 = wx_hi*wx_hi + wy_lo*wy_lo;
                    const Real djk_10 = wx_lo*wx_lo + wy_hi*wy_hi;
                    const Real djk_11 = wx_lo*wx_lo + wy_lo*wy_lo;

                    const Real wjk_00 = mij_00 * djk_01 * djk_10 * djk_11;
                    const Real wjk_01 = djk_00 * mij_01 * djk_10 * djk_11;
                    const Real wjk_10 = djk_00 * djk_01 * mij_10 * djk_11;
                    const Real wjk_11 = djk_00 * djk_01 * djk_10 * mij_11;

                    const Real inv_denom = 1.0 / (wjk_00 + wjk_01 + wjk_10 + wjk_11);

                    ep_s_loc += 3.0*(
                      + wjk_00*eps_arr[0][0]
                      + wjk_01*eps_arr[0][1]
                      + wjk_10*eps_arr[1][0]
                      + wjk_11*eps_arr[1][1])*inv_denom;

                    dtaudz = dxi[2]*(
                      + wjk_00*dtau_arr[0][0]
                      + wjk_01*dtau_arr[0][1]
                      + wjk_10*dtau_arr[1][0]
                      + wjk_11*dtau_arr[1][1])*inv_denom;

                  } else if ( count == 2 ) {

                    ep_s_loc += 2.0 * (
                      + mij_00*mij_10*(wx_lo*eps_arr[0][0] + wx_hi*eps_arr[1][0])
                      + mij_00*mij_01*(wy_lo*eps_arr[0][0] + wy_hi*eps_arr[0][1])
                      + mij_11*mij_10*(wy_lo*eps_arr[1][0] + wy_hi*eps_arr[1][1])
                      + mij_11*mij_01*(wx_lo*eps_arr[0][1] + wx_hi*eps_arr[1][1]));

                    dtaudz = dxi[2]*(
                      + mij_00*mij_10*(wx_lo*dtau_arr[0][0] + wx_hi*dtau_arr[1][0])
                      + mij_00*mij_01*(wy_lo*dtau_arr[0][0] + wy_hi*dtau_arr[0][1])
                      + mij_11*mij_10*(wy_lo*dtau_arr[1][0] + wy_hi*dtau_arr[1][1])
                      + mij_11*mij_01*(wx_lo*dtau_arr[0][1] + wx_hi*dtau_arr[1][1]));

                  } else {

                    AMREX_ALWAYS_ASSERT( count == 1);

                    ep_s_loc +=
                      + mij_00*eps_arr[0][0]
                      + mij_01*eps_arr[0][1]
                      + mij_10*eps_arr[1][0]
                      + mij_11*eps_arr[1][1];

                    dtaudz = dxi[2]*(
                      + mij_00*dtau_arr[0][0]
                      + mij_01*dtau_arr[0][1]
                      + mij_10*dtau_arr[1][0]
                      + mij_11*dtau_arr[1][1]);
                  }

                } // END dtau_dz

                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_wt > 0.,"Invalid ep_s_wt");
                ep_s_loc /= ep_s_wt;

              } // Cut cell


              AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_loc >= 0.,"ep_s_loc too low");
              AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_loc <= 1.,"ep_s_loc too high");

              // Store local solids volume fraction in unused array location.
              ptile_data.m_runtime_rdata[runtime_idxs.ep_s][pid] = ep_s_loc;


              p_realarray[SoArealData::omegax][pid] = dtaudx;
              p_realarray[SoArealData::omegay][pid] = dtaudy;
              p_realarray[SoArealData::omegaz][pid] = dtaudz;
            }

          }); // part loop
        } // if box not all regular

        /********************************************************************
         * Update runtime cost (used in load-balancing)                     *
         *******************************************************************/
        if (a_loadbalance->WeightByParticleRunTime()) {

          Real const weight( ParallelDescriptor::second() - timer_start );
          a_loadbalance->UpdateWeights(lev, pti.index(), weight);

        } else if (a_loadbalance->WeightByParticleCount()) {

          Real const weight(pti.numParticles() );
          a_loadbalance->UpdateWeights(lev, pti.index(), weight);

        }

      } // pti
    }  // omp region
  }// lev
}
