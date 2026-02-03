#ifndef _MFIX_EB_INTERP_SHEPARD_K_H_
#define _MFIX_EB_INTERP_SHEPARD_K_H_

#include <mfix_pc.H>

/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void shepard_interp (const amrex::RealVect& pos,
                     const int i,
                     const int j,
                     const int k,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
                     amrex::Array4<const amrex::EBCellFlag> const& flags,
                     amrex::Array4<const amrex::Real> const& ccent,
                     amrex::Array4<const amrex::Real> const& bcent,
                     amrex::Array4<const amrex::Real> const& apx,
                     amrex::Array4<const amrex::Real> const& apy,
                     amrex::Array4<const amrex::Real> const& apz,
                     amrex::Array4<amrex::Real const> const& interp_array,
                     amrex::Real* interp_vals,
                     const int srccomp,
                     const int dstcomp,
                     const int numcomp)
{

    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

    const amrex::Real cdist_x = pos[0] - (i + 0.5 + ccent(i,j,k,0))*dx[0] - plo[0];
    const amrex::Real cdist_y = pos[1] - (j + 0.5 + ccent(i,j,k,1))*dx[1] - plo[1];
    const amrex::Real cdist_z = pos[2] - (k + 0.5 + ccent(i,j,k,2))*dx[2] - plo[2];

    // Distance between particle and cell centoid.
    const amrex::Real cdist =
      std::sqrt(cdist_x*cdist_x + cdist_y*cdist_y + cdist_z*cdist_z);

    // Before doing anything fance, just check that the particle isn't overlapping
    // the cell center. If it is, use the cell value.
    if(cdist < tolerance) {

        for(int n=0; n<numcomp; n++)
            interp_vals[dstcomp+n] = interp_array(i,j,k,srccomp+n);

        return;
    }

    amrex::Real par_dot_EB(2.);

    // If the particle is inside a cut-cell, verify that it is on the correct
    // side of the EB before trying to interpolate.
    if(flags(i,j,k).isSingleValued()){

        const amrex::Real adx = apx(i+1,j  ,k  ) - apx(i  ,j  ,k  );
        const amrex::Real ady = apy(i  ,j+1,k  ) - apy(i  ,j  ,k  );
        const amrex::Real adz = apz(i  ,j  ,k+1) - apz(i  ,j  ,k  );

        const amrex::Real apnorm = std::sqrt(adx*adx + ady*ady + adz*adz);

        // Temporary sanity check.
        AMREX_ASSERT_WITH_MESSAGE( apnorm > std::numeric_limits<amrex::Real>::epsilon(),
            "Invalid EB normal. Not sure how to interpolate!");

        // EB normal
        const amrex::Real apnorminv = 1. / apnorm;
        const amrex::RealVect normal = {adx*apnorminv, ady*apnorminv, adz*apnorminv};

        // Projection of vector pointing from EB centroid to particle onto EB normal
        par_dot_EB = (pos[0] - (i + 0.5 + bcent(i,j,k,0))*dx[0] - plo[0])*normal[0] +
                     (pos[1] - (j + 0.5 + bcent(i,j,k,1))*dx[1] - plo[1])*normal[1] +
                     (pos[2] - (k + 0.5 + bcent(i,j,k,2))*dx[2] - plo[2])*normal[2];

#ifdef AMREX_DEBUG
        amrex::Real cent_dot_EB(1.);
        // Projection of vector pointing from EB centroid to cell centroid onto EB normal
        cent_dot_EB = (ccent(i,j,k,0) - bcent(i,j,k,0))*dx[0]*normal[0] +
                      (ccent(i,j,k,1) - bcent(i,j,k,1))*dx[1]*normal[1] +
                      (ccent(i,j,k,2) - bcent(i,j,k,2))*dx[2]*normal[2];

        // Temporary sanity check
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( cent_dot_EB > tolerance,
            "cent_dot_EB < tolerance ... this makes no sense!");
#endif
    }

    // A negative value implies that the particle is 'behind' the EB (and
    // therefore inside the wall). Although this shouldn't occur often, it could
    // for fast moving particles.  Set the fluid velocity to zero and use the
    // cell value for remaining items.
    if( flags(i,j,k).isSingleValued() && (par_dot_EB <= tolerance)){

      for(int n=0; n<numcomp; n++)
        interp_vals[dstcomp+n] = interp_array(i,j,k,srccomp+n);

    // The particle is near the EB. It is either
    // 1) between the cell centroid and the EB, or
    // 2) one or more cells in the 2^3 stencil is covered.
    // Either way, interpolating to the particle is not straight forward so
    // (for now) we do a 1D interpolation to the EB.

    } else {

      // The particle is in cell (i,j,k). Compute the non dimensional
      // distance from the cell center.

      const amrex::RealVect pcent = {
         (pos[0] - plo[0])*dxi[0] - static_cast<amrex::Real>(i) - 0.5,
         (pos[1] - plo[1])*dxi[1] - static_cast<amrex::Real>(j) - 0.5,
         (pos[2] - plo[2])*dxi[2] - static_cast<amrex::Real>(k) - 0.5};

      constexpr int wt_size = 27;

      amrex::Array1D<amrex::Real,0,wt_size-1> wt;
      amrex::Real max_wt(-1.);

      for(int kk = k-1; kk <= k+1; kk++) {
        for(int jj = j-1; jj <= j+1; jj++) {
          for(int ii = i-1; ii <= i+1; ii++) {

            const int cc_ind = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));

            AMREX_ALWAYS_ASSERT( 0 <= cc_ind && cc_ind < wt_size);

            wt(cc_ind) = 100.;

            if(!flags(ii,jj,kk).isCovered()) {

              amrex::Real x_off = static_cast<amrex::Real>(ii-i) - pcent[0];
              amrex::Real y_off = static_cast<amrex::Real>(jj-j) - pcent[1];
              amrex::Real z_off = static_cast<amrex::Real>(kk-k) - pcent[2];

              { // Distance from particle to cell centroid

                amrex::Real lx = x_off + ccent(ii,jj,kk,0);
                amrex::Real ly = y_off + ccent(ii,jj,kk,1);
                amrex::Real lz = z_off + ccent(ii,jj,kk,2);

                wt(cc_ind) = std::sqrt(lx*lx + ly*ly + lz*lz);
                max_wt = amrex::max(max_wt, wt(cc_ind));

              }

            } // end not covered cell
          } // end kk loop
        } // end jj loop
      } // end ii loop

      amrex::Real sum_wt(0.);
      for(int lc=0; lc<wt_size; lc++) {
        if (wt(lc) < std::numeric_limits<double>::min() ) {
          wt(lc) = 1.0e8;
          sum_wt += wt(lc);
        } else if (wt(lc) < max_wt) {
          wt(lc) = (max_wt - wt(lc))/(max_wt*wt(lc));
          wt(lc) *= wt(lc);
          sum_wt += wt(lc);
        } else {
          wt(lc) = 0.0;
        }
      }

      AMREX_ALWAYS_ASSERT(sum_wt > 0.0);

      const amrex::Real inv_sum_wt = 1.0/sum_wt;
      for(int lc=0; lc<wt_size; lc++) {
        wt(lc) *= inv_sum_wt;
      }


      for(int nn = 0; nn < numcomp; nn++ ) {

        // Clear array for accumulation.
        interp_vals[dstcomp+nn] = 0.0;

        for(int kk = k-1; kk <= k+1; kk++) {
          for(int jj = j-1; jj <= j+1; jj++) {
            for(int ii = i-1; ii <= i+1; ii++) {

              const int cc_ind = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));

              AMREX_ALWAYS_ASSERT( 0 <= cc_ind && cc_ind < wt_size);

              if(!flags(ii,jj,kk).isCovered())
                interp_vals[dstcomp+nn] += wt(cc_ind) * interp_array(ii,jj,kk,srccomp+nn);

            } // end ii loop
          } // end jj loop
        } // end kk loop
      } // end nn loop
    }
}


/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void shepard_interp_eb(const amrex::RealVect& pos,
                       const int i,
                       const int j,
                       const int k,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
                       amrex::Array4<const amrex::EBCellFlag> const& flags,
                       amrex::Array4<const amrex::Real> const& ccent,
                       amrex::Array4<const amrex::Real> const& bcent,
                       amrex::Array4<const amrex::Real> const& apx,
                       amrex::Array4<const amrex::Real> const& apy,
                       amrex::Array4<const amrex::Real> const& apz,
                       amrex::Array4<amrex::Real const> const& interp_array,
                       amrex::Real* interp_vals,
                       const int srccomp,
                       const int dstcomp,
                       const int numcomp)
{

    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

    const amrex::Real cdist_x = pos[0] - (i + 0.5 + ccent(i,j,k,0))*dx[0] - plo[0];
    const amrex::Real cdist_y = pos[1] - (j + 0.5 + ccent(i,j,k,1))*dx[1] - plo[1];
    const amrex::Real cdist_z = pos[2] - (k + 0.5 + ccent(i,j,k,2))*dx[2] - plo[2];

    // Distance between particle and cell centoid.
    const amrex::Real cdist =
      std::sqrt(cdist_x*cdist_x + cdist_y*cdist_y + cdist_z*cdist_z);

    // Before doing anything fance, just check that the particle isn't overlapping
    // the cell center. If it is, use the cell value.
    if(cdist < tolerance) {

        for(int n=0; n<numcomp; n++)
            interp_vals[dstcomp+n] = interp_array(i,j,k,srccomp+n);

        return;
    }

    amrex::Real par_dot_EB(2.);

    // If the particle is inside a cut-cell, verify that it is on the correct
    // side of the EB before trying to interpolate.
    if(flags(i,j,k).isSingleValued()){

        const amrex::Real adx = apx(i+1,j  ,k  ) - apx(i  ,j  ,k  );
        const amrex::Real ady = apy(i  ,j+1,k  ) - apy(i  ,j  ,k  );
        const amrex::Real adz = apz(i  ,j  ,k+1) - apz(i  ,j  ,k  );

        const amrex::Real apnorm = std::sqrt(adx*adx + ady*ady + adz*adz);

        // Temporary sanity check.
        AMREX_ASSERT_WITH_MESSAGE( apnorm > std::numeric_limits<amrex::Real>::epsilon(),
            "Invalid EB normal. Not sure how to interpolate!");

        // EB normal
        const amrex::Real apnorminv = 1. / apnorm;
        const amrex::RealVect normal = {adx*apnorminv, ady*apnorminv, adz*apnorminv};

        // Projection of vector pointing from EB centroid to particle onto EB normal
        par_dot_EB = (pos[0] - (i + 0.5 + bcent(i,j,k,0))*dx[0] - plo[0])*normal[0] +
                     (pos[1] - (j + 0.5 + bcent(i,j,k,1))*dx[1] - plo[1])*normal[1] +
                     (pos[2] - (k + 0.5 + bcent(i,j,k,2))*dx[2] - plo[2])*normal[2];

#ifdef AMREX_DEBUG

        amrex::Real cent_dot_EB(1.);
        // Projection of vector pointing from EB centroid to cell centroid onto EB normal
        cent_dot_EB = (ccent(i,j,k,0) - bcent(i,j,k,0))*dx[0]*normal[0] +
          (ccent(i,j,k,1) - bcent(i,j,k,1))*dx[1]*normal[1] +
          (ccent(i,j,k,2) - bcent(i,j,k,2))*dx[2]*normal[2];

        // Temporary sanity check
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( cent_dot_EB > tolerance,
            "cent_dot_EB < tolerance ... this makes no sense!");
#endif


    }

    // A negative value implies that the particle is 'behind' the EB (and
    // therefore inside the wall). Although this shouldn't occur often, it could
    // for fast moving particles.  Set the fluid velocity to zero and use the
    // cell value for remaining items.
    if( flags(i,j,k).isSingleValued() && (par_dot_EB <= tolerance)){

      for(int n=0; n<numcomp; n++)
        interp_vals[dstcomp+n] = interp_array(i,j,k,srccomp+n);

    // The particle is near the EB. It is either
    // 1) between the cell centroid and the EB, or
    // 2) one or more cells in the 2^3 stencil is covered.
    // Either way, interpolating to the particle is not straight forward so
    // (for now) we do a 1D interpolation to the EB.

    } else {

      // The particle is in cell (i,j,k). Compute the non dimensional
      // distance from the cell center.

      const amrex::RealVect pcent = {
         (pos[0] - plo[0])*dxi[0] - static_cast<amrex::Real>(i) - 0.5,
         (pos[1] - plo[1])*dxi[1] - static_cast<amrex::Real>(j) - 0.5,
         (pos[2] - plo[2])*dxi[2] - static_cast<amrex::Real>(k) - 0.5};

      constexpr int wt_size = 54;
      amrex::Array1D<amrex::Real,0,wt_size-1> wt;
      amrex::Real max_wt(-1.);

      for(int kk = k-1; kk <= k+1; kk++) {
        for(int jj = j-1; jj <= j+1; jj++) {
          for(int ii = i-1; ii <= i+1; ii++) {

            const int cc_ind = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));
            const int eb_ind = cc_ind + 27;

            AMREX_ALWAYS_ASSERT( 0 <= cc_ind && cc_ind <= 26);
            AMREX_ALWAYS_ASSERT(27 <= eb_ind && eb_ind <= 53);

            wt(cc_ind) = 100.;
            wt(eb_ind) = 100.;

            if(!flags(ii,jj,kk).isCovered()) {

              amrex::Real x_off = static_cast<amrex::Real>(ii-i) - pcent[0];
              amrex::Real y_off = static_cast<amrex::Real>(jj-j) - pcent[1];
              amrex::Real z_off = static_cast<amrex::Real>(kk-k) - pcent[2];

              { // Distance from particle to cell centroid

                amrex::Real lx = x_off + ccent(ii,jj,kk,0);
                amrex::Real ly = y_off + ccent(ii,jj,kk,1);
                amrex::Real lz = z_off + ccent(ii,jj,kk,2);

                wt(cc_ind) = std::sqrt(lx*lx + ly*ly + lz*lz);
                max_wt = amrex::max(max_wt, wt(cc_ind));

              }

              if (!flags(ii,jj,kk).isRegular()) {
                // Distance from particle to eb centroid

                amrex::Real lx = x_off + bcent(ii,jj,kk,0);
                amrex::Real ly = y_off + bcent(ii,jj,kk,1);
                amrex::Real lz = z_off + bcent(ii,jj,kk,2);

                wt(eb_ind) = std::sqrt(lx*lx + ly*ly + lz*lz);
                max_wt = amrex::max(max_wt, wt(eb_ind));
              }
            } // end not covered cell
          } // end kk loop
        } // end jj loop
      } // end ii loop

      amrex::Real sum_wt(0.);
      for(int lc=0; lc<wt_size; lc++) {
        if (wt(lc) < std::numeric_limits<double>::min() ) {
          wt(lc) = 1.0e8;
          sum_wt += wt(lc);
        } else if (wt(lc) < max_wt) {
          wt(lc) = (max_wt - wt(lc))/(max_wt*wt(lc));
          wt(lc) *= wt(lc);
          sum_wt += wt(lc);
        } else {
          wt(lc) = 0.0;
        }
      }

      AMREX_ALWAYS_ASSERT(sum_wt > 0.0);

      const amrex::Real inv_sum_wt = 1.0/sum_wt;
      for(int lc=0; lc<wt_size; lc++) {
        wt(lc) *= inv_sum_wt;
      }


      for(int nn = 0; nn < numcomp; nn++ ) {

        // Clear array for accumulation.
        interp_vals[dstcomp+nn] = 0.0;

        for(int kk = k-1; kk <= k+1; kk++) {
          for(int jj = j-1; jj <= j+1; jj++) {
            for(int ii = i-1; ii <= i+1; ii++) {

              const int cc_ind = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));
              const int eb_ind = cc_ind + 27;

              AMREX_ALWAYS_ASSERT( 0 <= cc_ind && cc_ind <= 26);
              AMREX_ALWAYS_ASSERT(27 <= eb_ind && eb_ind <= 53);

              if(!flags(ii,jj,kk).isCovered()) {

                interp_vals[dstcomp+nn] += wt(cc_ind) * interp_array(ii,jj,kk,srccomp+nn);

                // Add in the value from the the EB.
                // HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK
                // Currently assuming that the EB has a value of ZERO for no slip walls.
                //if (!flags(ii,jj,kk).isRegular()) {
                //  interp_vals[dstcomp+n] += wt(eb_ind) * --EB Diriclet value--;
                //}
              }
            } // end ii loop
          } // end jj loop
        } // end kk loop
      } // end nn loop
    }
}

#endif
