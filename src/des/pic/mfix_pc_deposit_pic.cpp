#include <AMReX.H>
#include "AMReX_Particles.H"
#include <mfix_pc.H>

#include <mfix_deposition_K.H>
#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_pic_K.H>

#include <mfix_bc.H>
#include <mfix_pic.H>

using namespace amrex;

void MFIXParticleContainer::
MFIX_PC_SolidsVelocityDeposition (int lev,
                                  Array<MultiFab*,3>& vel_s_in,
                                  const FabArray<EBCellFlagFab>* flags)
{
  BL_PROFILE("(MFIXParticleContainer::MFIX_PC_SolidsVelocityDeposition)");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto p_lo = gm.ProbLoArray();
  const auto  dxi = gm.InvCellSizeArray();

  const auto inv_reg_cell_vol = dxi[0]*dxi[1]*dxi[2];

  {
    FArrayBox local_u_s_fab;
    FArrayBox local_v_s_fab;
    FArrayBox local_w_s_fab;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();

      FArrayBox& u_s_fab = (*vel_s_in[0])[pti];
      FArrayBox& v_s_fab = (*vel_s_in[1])[pti];
      FArrayBox& w_s_fab = (*vel_s_in[2])[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered ) {

        auto u_s = u_s_fab.array();
        auto v_s = v_s_fab.array();
        auto w_s = w_s_fab.array();


        amrex::ParallelFor(nrp, [pstruct,p_realarray,p_lo,dxi,inv_reg_cell_vol,u_s, v_s, w_s]
        AMREX_GPU_DEVICE (int ip) noexcept
        {
          const ParticleType& p = pstruct[ip];

          const Real lx = (p.pos(0) - p_lo[0]) * dxi[0] + 0.5;
          const Real ly = (p.pos(1) - p_lo[1]) * dxi[1] + 0.5;
          const Real lz = (p.pos(2) - p_lo[2]) * dxi[2] + 0.5;

          const int i = static_cast<int>(amrex::Math::floor(lx));
          const int j = static_cast<int>(amrex::Math::floor(ly));
          const int k = static_cast<int>(amrex::Math::floor(lz));

          const Real wx_hi(lx - static_cast<Real>(i));
          const Real wy_hi(ly - static_cast<Real>(j));
          const Real wz_hi(lz - static_cast<Real>(k));

          const Real wx_lo(1.0 - wx_hi);
          const Real wy_lo(1.0 - wy_hi);
          const Real wz_lo(1.0 - wz_hi);

          const Real pvol = p_realarray[SoArealData::statwt][ip] *
            SoArealData::volume(p_realarray[SoArealData::radius][ip]) *
            inv_reg_cell_vol;

          {// Deposition of x velocity -- x-face deposition

            const Real pvelx = pvol*p_realarray[SoArealData::velx][ip];

            HostDevice::Atomic::Add(&u_s(i   , j-1, k-1, 0),wy_lo*wz_lo*pvelx);
            HostDevice::Atomic::Add(&u_s(i   , j-1, k  , 0),wy_lo*wz_hi*pvelx);
            HostDevice::Atomic::Add(&u_s(i   , j,   k-1, 0),wy_hi*wz_lo*pvelx);
            HostDevice::Atomic::Add(&u_s(i   , j,   k  , 0),wy_hi*wz_hi*pvelx);

            HostDevice::Atomic::Add(&u_s(i   , j-1, k-1, 1),wy_lo*wz_lo*pvol);
            HostDevice::Atomic::Add(&u_s(i   , j-1, k  , 1),wy_lo*wz_hi*pvol);
            HostDevice::Atomic::Add(&u_s(i   , j,   k-1, 1),wy_hi*wz_lo*pvol);
            HostDevice::Atomic::Add(&u_s(i   , j,   k  , 1),wy_hi*wz_hi*pvol);
          }


          {// Deposition of y velocity -- y-face deposition

            const Real pvely = pvol*p_realarray[SoArealData::vely][ip];

            HostDevice::Atomic::Add(&v_s(i-1, j ,   k-1, 0),wx_lo*wz_lo*pvely);
            HostDevice::Atomic::Add(&v_s(i-1, j ,   k  , 0),wx_lo*wz_hi*pvely);
            HostDevice::Atomic::Add(&v_s(i,   j ,   k-1, 0),wx_hi*wz_lo*pvely);
            HostDevice::Atomic::Add(&v_s(i,   j ,   k  , 0),wx_hi*wz_hi*pvely);

            HostDevice::Atomic::Add(&v_s(i-1, j ,   k-1, 1),wx_lo*wz_lo*pvol);
            HostDevice::Atomic::Add(&v_s(i-1, j ,   k  , 1),wx_lo*wz_hi*pvol);
            HostDevice::Atomic::Add(&v_s(i,   j ,   k-1, 1),wx_hi*wz_lo*pvol);
            HostDevice::Atomic::Add(&v_s(i,   j ,   k  , 1),wx_hi*wz_hi*pvol);
          }


          {// Deposition of z velocity -- z-face deposition

            const Real pvelz = pvol*p_realarray[SoArealData::velz][ip];

            HostDevice::Atomic::Add(&w_s(i-1, j-1, k   , 0),wx_lo*wy_lo*pvelz);
            HostDevice::Atomic::Add(&w_s(i-1, j,   k   , 0),wx_lo*wy_hi*pvelz);
            HostDevice::Atomic::Add(&w_s(i,   j-1, k   , 0),wx_hi*wy_lo*pvelz);
            HostDevice::Atomic::Add(&w_s(i,   j,   k   , 0),wx_hi*wy_hi*pvelz);

            HostDevice::Atomic::Add(&w_s(i-1, j-1, k   , 1),wx_lo*wy_lo*pvol);
            HostDevice::Atomic::Add(&w_s(i-1, j,   k   , 1),wx_lo*wy_hi*pvol);
            HostDevice::Atomic::Add(&w_s(i,   j-1, k   , 1),wx_hi*wy_lo*pvol);
            HostDevice::Atomic::Add(&w_s(i,   j,   k   , 1),wx_hi*wy_hi*pvol);
          }
        });

      }
    }
  }
}


void
MFIXParticleContainer::PICHydroStep (int lev,
                                     const bool apply_forces,
                                     const bool update_parcels,
                                     const bool use_taylor_approx,
                                     const Real advance_vel_p,
                                     Real dt,
                                     RealVect& gravity,
                                     Vector< Array<MultiFab*,3> >& vel_s_in,
                                     MultiFab & ep_s_out,
                                     Array<MultiFab*,3>& vel_s_out,
                                     const MultiFab * volfrac,
                                     const amrex::FabArray<EBCellFlagFab>* flags,
                                     EBFArrayBoxFactory* ebfactory,
                                     const int ls_refinement,
                                     const MultiFab* ls_phi)
{
  BL_PROFILE("MFIXParticleContainer::SolidsVolumeDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      p_lo = gm.ProbLoArray();
  const auto      p_hi = gm.ProbHiArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  GpuArray<int,3> periodic = gm.isPeriodicArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

  Real const en = (m_pic.damping_factor() + 1.0);

  Real const en_w = m_pic.damping_factor_wall_normal();
  Real const et_w = m_pic.damping_factor_wall_tangent();

  Real const vel_ref_frame = m_pic.vel_ref_frame();

  Real const ep_cp = m_pic.ep_cp();
  Real const inv_ep_cp = 1.0/ep_cp;

  Real const three_sqrt_two(3.0*std::sqrt(2.0));

  int const x_lo_bc = m_boundary_conditions.domain_bc(0);
  int const x_hi_bc = m_boundary_conditions.domain_bc(1);
  int const y_lo_bc = m_boundary_conditions.domain_bc(2);
  int const y_hi_bc = m_boundary_conditions.domain_bc(3);
  int const z_lo_bc = m_boundary_conditions.domain_bc(4);
  int const z_hi_bc = m_boundary_conditions.domain_bc(5);

  auto& plev  = GetParticles(lev);

  int const include_vm( m_runtimeRealData.contains_vm() );
  int const idx_pc_vm_coeff( m_runtimeRealData.vm_coeff );

  if (include_vm) {
    AMREX_ASSERT( idx_pc_vm_coeff >= 0);
    AMREX_ASSERT( idx_pc_vm_coeff < m_runtimeRealData.count );
  }

  int const include_acc( m_runtimeRealData.contains_acc() );
  int const idx_pc_acc( m_runtimeRealData.acceleration );
  if (include_acc) {
    AMREX_ASSERT( idx_pc_acc >= 0);
    AMREX_ASSERT( idx_pc_acc < m_runtimeRealData.count );
  }

  {
    FArrayBox local_fab_to_be_filled;

    FArrayBox local_u_s_fab;
    FArrayBox local_v_s_fab;
    FArrayBox local_w_s_fab;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& ptile = plev[index];
      auto ptile_data = ptile.getParticleTileData();
      auto runtime_idxs = m_runtimeRealData;

      auto& particles = pti.GetArrayOfStructs();
      ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();
      FArrayBox& fab_to_be_filled = ep_s_out[pti];

      FArrayBox& u_s_fab = (*vel_s_out[0])[pti];
      FArrayBox& v_s_fab = (*vel_s_out[1])[pti];
      FArrayBox& w_s_fab = (*vel_s_out[2])[pti];

      const Box& bx  = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(bx) != FabType::covered ) {

        auto volarr = fab_to_be_filled.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto& vfrac = (*volfrac)[pti].array();

        // Determine if this particle tile actually has any walls
        bool has_walls = false;

        if ((ebfactory != NULL) &&
           ((*flags)[pti].getType(amrex::grow(bx,1)) == FabType::singlevalued))  {

          has_walls = true;

        } else {
          // We need this test for the case of an inflow boundary:
          // inflow does not appear in the EBFactory but
          // the particles see it as a wall

          // Create the nodal refined box based on the current particle tile
          Box refined_box(amrex::convert(amrex::refine(bx,ls_refinement), IntVect{1,1,1}));

          // Set tol to 1/2 dx
          Real tol = amrex::min(dx[0], amrex::min(dx[1], dx[2])) / 2;

          Real ls_min_over_box = ((*ls_phi)[pti]).min<RunOn::Gpu>(refined_box,0);

          if (ls_min_over_box < tol) has_walls = true;

        }

        const auto& phiarr = (has_walls) ? ls_phi->const_array(pti) : Array4<Real const>{};

        //const auto& avg_prop_array = avg_prop_in[lev]->array(pti);
        const auto& u_so = (*vel_s_in[lev][0]).const_array(pti);
        const auto& v_so = (*vel_s_in[lev][1]).const_array(pti);
        const auto& w_so = (*vel_s_in[lev][2]).const_array(pti);

        auto u_s = u_s_fab.array();
        auto v_s = v_s_fab.array();
        auto w_s = w_s_fab.array();

        amrex::ParallelFor(nrp,
           [pstruct,p_realarray,p_hi,p_lo,dx,dxi,vfrac,volarr, u_so, v_so, w_so, en, ep_cp,
            reg_cell_vol,flagsarr, dt, gravity, has_walls,ls_refinement,phiarr,
            vel_ref_frame, three_sqrt_two, en_w, et_w, u_s, v_s, w_s, inv_ep_cp,
            apply_forces, update_parcels, use_taylor_approx, advance_vel_p,
            x_lo_bc,x_hi_bc, y_lo_bc,y_hi_bc,z_lo_bc,z_hi_bc, periodic, ptile_data,
            include_vm, idx_pc_vm_coeff, include_acc, idx_pc_acc, runtime_idxs]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            ParticleType& p = pstruct[ip];

            const amrex::RealVect vel_p_old  = {p_realarray[SoArealData::velx][ip],
                                                p_realarray[SoArealData::vely][ip],
                                                p_realarray[SoArealData::velz][ip]};

            // position
            RealVect pos(p.pos());
            RealVect vel_p(vel_p_old);
            Real const vol_p(SoArealData::volume(p_realarray[SoArealData::radius][ip]));

            if (apply_forces ) {

              amrex::Real const mass_p = p_realarray[SoArealData::density][ip]*vol_p;

              // coeff multiplies 'old velocity' in update
              amrex::Real coeff(mass_p);
              if (include_vm) { coeff += ptile_data.m_runtime_rdata[idx_pc_vm_coeff][ip]; }

              // update denominator
              Real const scale( 1.0 / (coeff + dt*p_realarray[SoArealData::drag_coeff][ip]) );

              // Sources acting on parcle: gravity, grad_p, div(tau), beta*u_f, Cvm*DufDt, ...
              const amrex::RealVect S_c = { mass_p*gravity[0] + p_realarray[SoArealData::vel_source_x][ip],
                                            mass_p*gravity[1] + p_realarray[SoArealData::vel_source_y][ip],
                                            mass_p*gravity[2] + p_realarray[SoArealData::vel_source_z][ip]};

              // solids volume fraction
              Real const eps_p( amrex::max(1.0e-8, ptile_data.m_runtime_rdata[runtime_idxs.ep_s][ip]) );

              // solids stress gradient:
              const amrex::RealVect grad_tau_p = {p_realarray[SoArealData::omegax][ip],
                                                  p_realarray[SoArealData::omegay][ip],
                                                  p_realarray[SoArealData::omegaz][ip]};

              const Real mfp_vel = (p_realarray[SoArealData::radius][ip] /
                                    (dt * three_sqrt_two * eps_p));

              // Compute the updated PIC velocity
              vel_p = updated_pic_velocity(pos, vel_p_old, dt, coeff, scale, S_c, grad_tau_p,
                 u_so, v_so, w_so, en, vol_p, eps_p, ep_cp, vel_ref_frame, mfp_vel, dxi, p_lo);

              // Update parcel positions
              pos[0] += dt * vel_p[0];
              pos[1] += dt * vel_p[1];
              pos[2] += dt * vel_p[2];

              // Effective radius of the parcel
              Real eff_radius = p_realarray[SoArealData::radius][ip] *
                std::cbrt(p_realarray[SoArealData::statwt][ip] * inv_ep_cp);

              // If this FAB has EB, reflect any parcels overlapping the levelset
              // back into the domain.
              if (has_walls) {

                Real ls_value = interp_level_set(pos, ls_refinement, phiarr, p_lo, dxi);

                const Real overlap = eff_radius - ls_value;

                // The particle intersects the levelset.
                if (overlap > 0.) {

                  RealVect normal(0.);
                  level_set_normal(pos, ls_refinement, normal, phiarr, p_lo, dxi);

                  // Reflect the particle.
                  pos[0] += overlap*normal[0];
                  pos[1] += overlap*normal[1];
                  pos[2] += overlap*normal[2];

                  // Plane ref point
                  const Real Nw_Vp = normal[0]*vel_p[0] + normal[1]*vel_p[1] + normal[2]*vel_p[2];

                  // Parcel normal velocity
                  const RealVect Vpn = {Nw_Vp*normal[0], Nw_Vp*normal[1], Nw_Vp*normal[2]};

                  // Parcel tangential velocity
                  const RealVect Vpt = {vel_p[0]-Vpn[0], vel_p[1]-Vpn[1], vel_p[2]-Vpn[2]};

                  // Rebound parcel if moving towards wall.
                  if(Nw_Vp < 0.) {
                    vel_p[0] = -en_w*Vpn[0] + et_w*Vpt[0];
                    vel_p[1] = -en_w*Vpn[1] + et_w*Vpt[1];
                    vel_p[2] = -en_w*Vpn[2] + et_w*Vpt[2];

                  } else {
                    vel_p[0] = Vpn[0] + et_w*Vpt[0];
                    vel_p[1] = Vpn[1] + et_w*Vpt[1];
                    vel_p[2] = Vpn[2] + et_w*Vpt[2];
                  }

                }
              } // end has_walls


              // Impose domain constraints.
              if (pos[0] < p_lo[0]+eff_radius) {
                if (x_lo_bc) { pos[0] = p_lo[0] + eff_radius; }
                else if (update_parcels && !periodic[0]) { p.id().make_invalid(); }
              }

              if (pos[0]+eff_radius > p_hi[0]) {
                if (x_hi_bc) { pos[0] = p_hi[0] - eff_radius; }
                else if (update_parcels && !periodic[0]) { p.id().make_invalid(); }
              }

              if (pos[1] < p_lo[1]+eff_radius) {
                if (y_lo_bc) { pos[1] = p_lo[1] + eff_radius; }
                else if (update_parcels && !periodic[1]) { p.id().make_invalid(); }
              }

              if (pos[1] + eff_radius > p_hi[1]) {
                if (y_hi_bc) { pos[1] = p_hi[1] - eff_radius; }
                else if (update_parcels && !periodic[1]) { p.id().make_invalid(); }
              }

              if (pos[2] < p_lo[2]+eff_radius) {
                if (z_lo_bc) { pos[2] = p_lo[2] + eff_radius; }
                else if (update_parcels && !periodic[2]) { p.id().make_invalid(); }
              }

              if (pos[2] + eff_radius > p_hi[2]) {
                if (z_hi_bc) { pos[2] = p_hi[2] - eff_radius; }
                else if (update_parcels && !periodic[2]) { p.id().make_invalid(); }
              }
            }

            if (p.id().is_valid()) {
              amrex::Real x = (pos[0] - p_lo[0]) * dxi[0] + 0.5;
              amrex::Real y = (pos[1] - p_lo[1]) * dxi[1] + 0.5;
              amrex::Real z = (pos[2] - p_lo[2]) * dxi[2] + 0.5;

              int i = static_cast<int>(amrex::Math::floor(x));
              int j = static_cast<int>(amrex::Math::floor(y));
              int k = static_cast<int>(amrex::Math::floor(z));

              amrex::GpuArray<amrex::Real,2> wx;
              amrex::GpuArray<amrex::Real,2> wy;
              amrex::GpuArray<amrex::Real,2> wz;

              wx[1] = x - static_cast<Real>(i);
              wy[1] = y - static_cast<Real>(j);
              wz[1] = z - static_cast<Real>(k);

              wx[0] = 1.0 - wx[1];
              wy[0] = 1.0 - wy[1];
              wz[0] = 1.0 - wz[1];

              amrex::Real total_weight = 0.0;
              GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

              if (use_taylor_approx) {

                const amrex::GpuArray<amrex::Real,2> ws({-1.0,1.0});

                const amrex::RealVect dt_vel_dxi =
                  {dt*p_realarray[SoArealData::velx][ip]*dxi[0],
                   dt*p_realarray[SoArealData::vely][ip]*dxi[1],
                   dt*p_realarray[SoArealData::velz][ip]*dxi[2]};

                for (int ii = 0; ii <= 1; ++ii)
                  for (int jj = 0; jj <= 1; ++jj)
                    for (int kk = 0; kk <= 1; ++kk){
                      if( !flagsarr(i-1+ii,j-1+jj,k-1+kk).isCovered() ) {

                        weights[ii][jj][kk] = wx[ii]*wy[jj]*wz[kk] +
                          ws[ii]*wy[jj]*wz[kk]*dt_vel_dxi[0] +
                          wx[ii]*ws[jj]*wz[kk]*dt_vel_dxi[1] +
                          wx[ii]*wy[jj]*ws[kk]*dt_vel_dxi[2];

                        total_weight += weights[ii][jj][kk];

                      } else {
                        weights[ii][jj][kk] = 0.0;
                      }
                    }

              } else {

                for (int ii = 0; ii <= 1; ++ii)
                  for (int jj = 0; jj <= 1; ++jj)
                    for (int kk = 0; kk <= 1; ++kk){
                      if( !flagsarr(i-1+ii,j-1+jj,k-1+kk).isCovered() ) {
                        weights[ii][jj][kk] = wx[ii]*wy[jj]*wz[kk];
                        total_weight += weights[ii][jj][kk];

                      } else {
                        weights[ii][jj][kk] = 0.0;
                      }
                    }
              }

              for (int ii = 0; ii <= 1; ++ii)
                for (int jj = 0; jj <= 1; ++jj)
                  for (int kk = 0; kk <= 1; ++kk)
                    weights[ii][jj][kk] /= total_weight;

              const Real pvol = p_realarray[SoArealData::statwt][ip] * vol_p / reg_cell_vol;

              for (int kk = -1; kk <= 0; ++kk) {
                for (int jj = -1; jj <= 0; ++jj) {
                  for (int ii = -1; ii <= 0; ++ii) {
                    if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                      continue;

                    Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);
                    HostDevice::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);
                  }
                }
              }

              const Real new_fac = advance_vel_p;
              const Real old_fac = 1.0 - new_fac;

              x = ((p.pos(0) + new_fac * dt * vel_p[0]) - p_lo[0]) * dxi[0] + 0.5;
              y = ((p.pos(1) + new_fac * dt * vel_p[1]) - p_lo[1]) * dxi[1] + 0.5;
              z = ((p.pos(2) + new_fac * dt * vel_p[2]) - p_lo[2]) * dxi[2] + 0.5;

              i = static_cast<int>(amrex::Math::floor(x));
              j = static_cast<int>(amrex::Math::floor(y));
              k = static_cast<int>(amrex::Math::floor(z));

              const Real wx_hi = x - static_cast<Real>(i);
              const Real wy_hi = y - static_cast<Real>(j);
              const Real wz_hi = z - static_cast<Real>(k);

              const Real wx_lo = 1.0 - wx_hi;
              const Real wy_lo = 1.0 - wy_hi;
              const Real wz_lo = 1.0 - wz_hi;


              {// Deposition of x velocity -- x-face deposition

                const Real pvelx = pvol*(new_fac*vel_p[0] + old_fac*vel_p_old[0]);

                HostDevice::Atomic::Add(&u_s(i   , j-1, k-1, 0),wy_lo*wz_lo*pvelx);
                HostDevice::Atomic::Add(&u_s(i   , j-1, k  , 0),wy_lo*wz_hi*pvelx);
                HostDevice::Atomic::Add(&u_s(i   , j,   k-1, 0),wy_hi*wz_lo*pvelx);
                HostDevice::Atomic::Add(&u_s(i   , j,   k  , 0),wy_hi*wz_hi*pvelx);

                HostDevice::Atomic::Add(&u_s(i   , j-1, k-1, 1),wy_lo*wz_lo*pvol);
                HostDevice::Atomic::Add(&u_s(i   , j-1, k  , 1),wy_lo*wz_hi*pvol);
                HostDevice::Atomic::Add(&u_s(i   , j,   k-1, 1),wy_hi*wz_lo*pvol);
                HostDevice::Atomic::Add(&u_s(i   , j,   k  , 1),wy_hi*wz_hi*pvol);
              }


              {// Deposition of y velocity -- y-face deposition

                const Real pvely = pvol*(new_fac*vel_p[1] + old_fac*vel_p_old[1]);

                HostDevice::Atomic::Add(&v_s(i-1, j ,   k-1, 0),wx_lo*wz_lo*pvely);
                HostDevice::Atomic::Add(&v_s(i-1, j ,   k  , 0),wx_lo*wz_hi*pvely);
                HostDevice::Atomic::Add(&v_s(i,   j ,   k-1, 0),wx_hi*wz_lo*pvely);
                HostDevice::Atomic::Add(&v_s(i,   j ,   k  , 0),wx_hi*wz_hi*pvely);

                HostDevice::Atomic::Add(&v_s(i-1, j ,   k-1, 1),wx_lo*wz_lo*pvol);
                HostDevice::Atomic::Add(&v_s(i-1, j ,   k  , 1),wx_lo*wz_hi*pvol);
                HostDevice::Atomic::Add(&v_s(i,   j ,   k-1, 1),wx_hi*wz_lo*pvol);
                HostDevice::Atomic::Add(&v_s(i,   j ,   k  , 1),wx_hi*wz_hi*pvol);
              }


              {// Deposition of z velocity -- z-face deposition

                const Real pvelz = pvol*(new_fac*vel_p[2] + old_fac*vel_p_old[2]);

                HostDevice::Atomic::Add(&w_s(i-1, j-1, k   , 0),wx_lo*wy_lo*pvelz);
                HostDevice::Atomic::Add(&w_s(i-1, j,   k   , 0),wx_lo*wy_hi*pvelz);
                HostDevice::Atomic::Add(&w_s(i,   j-1, k   , 0),wx_hi*wy_lo*pvelz);
                HostDevice::Atomic::Add(&w_s(i,   j,   k   , 0),wx_hi*wy_hi*pvelz);

                HostDevice::Atomic::Add(&w_s(i-1, j-1, k   , 1),wx_lo*wy_lo*pvol);
                HostDevice::Atomic::Add(&w_s(i-1, j,   k   , 1),wx_lo*wy_hi*pvol);
                HostDevice::Atomic::Add(&w_s(i,   j-1, k   , 1),wx_hi*wy_lo*pvol);
                HostDevice::Atomic::Add(&w_s(i,   j,   k   , 1),wx_hi*wy_hi*pvol);
              }


              if (update_parcels) {

                // Update positions
                p.pos(0) = pos[0];
                p.pos(1) = pos[1];
                p.pos(2) = pos[2];

                if (include_acc) {
                  ptile_data.m_runtime_rdata[idx_pc_acc  ][ip] =
                      (vel_p[0] - p_realarray[SoArealData::velx][ip])/dt;
                  ptile_data.m_runtime_rdata[idx_pc_acc+1][ip] =
                      (vel_p[1] - p_realarray[SoArealData::vely][ip])/dt;
                  ptile_data.m_runtime_rdata[idx_pc_acc+2][ip] =
                      (vel_p[2] - p_realarray[SoArealData::velz][ip])/dt;
                }

                // update parcel velocity
                p_realarray[SoArealData::velx][ip] = vel_p[0];
                p_realarray[SoArealData::vely][ip] = vel_p[1];
                p_realarray[SoArealData::velz][ip] = vel_p[2];
              }
            }
          });

        Gpu::synchronize();


      }
    }
  }
}
