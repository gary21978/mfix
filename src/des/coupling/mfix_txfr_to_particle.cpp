#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>

#include <mfix.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_shepard_K.H>
#include <mfix_filcc.H>
#include <mfix_deposition_op.H>
#include <mfix_mf_helpers.H>

void
mfix::
calc_txfr_particle ( Real const a_time, Real const a_dt,
                     Vector< MultiFab      *> const& a_vel,
                     Vector< MultiFab      *> const& a_T,
                     Vector< MultiFab const*> const& a_gradp,
                     Vector< MultiFab const*> const& a_divtau,
                     Vector< MultiFab const*> const& a_DufDt)
{
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::calc_txfr_particle()");

  //***************************************************************************
  //
  //***************************************************************************

  // Set flags for what to include in drag force.
  int const include_divtau( m_coupling.include_divtau() );
  int const include_vm( m_coupling.include_virtual_mass() );

  // Flag to include temperature
  int const include_Tf( fluid.solve_enthalpy() );

  GridToParticleIndexes::TxfrToParicle interp_idx( include_Tf,
      include_divtau, include_vm);

  for (int lev(0); lev < nlev(); lev++) {
    // Extrapolate velocity Dirichlet bc's to ghost cells
    bcs().set_velocity_bcs(lev, a_time, a_vel[lev], BCType::hoextrap);

    if (include_Tf) {
      bcs().set_temperature_bcs(lev, a_time, fluid, a_T[lev]);
    }
  }

  // This equals the biggest value interp_idx.count can return.
  // velocity +3, pressure gradient +3, temperature +1,
  // divtau +3; DufDt +3
  int constexpr interp_array_size(13);
  AMREX_ALWAYS_ASSERT( interp_idx.count <= interp_array_size );

  for (int lev(0); lev<nlev(); ++lev) {

    Box domain(geom[lev].Domain());

    MultiFab gp_tmp;
    gp_tmp.define(grids[lev], dmap[lev], 3, 1, MFInfo(), *(m_eb->factory()[lev]));

    MultiFab::Copy(gp_tmp, *a_gradp[lev], 0, 0, 3, 1);
    gp_tmp.FillBoundary(geom[lev].periodicity());

    for (MFIter mfi(gp_tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      // Call set_gradp_bcs after calling FillBoundary because the set_gradp_bcs
      // routine sets the ghost cells exterior to the domain from ghost cells
      // interior to the domain.
      m_boundary_conditions.set_gradp_bcs(bx, lev, gp_tmp[mfi], domain);
    }

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    // Pointer to Multifab for interpolation
    MultiFab* interp_ptr;

    // This is just a sanity check to make sure we're not using covered values
    // We can remove these lines once we're confident in the algorithm
    EB_set_covered(*a_vel[lev], 0, 3, 1, covered_val);
    EB_set_covered(gp_tmp, 0, 3, 1, covered_val);
    if (include_Tf) { EB_set_covered(*a_T[lev], 0, 1, 1, covered_val); }

    // Only one layer needed for interpolation
    const int interp_ng = 1;

    if (OnSameGrids) {

      // Store gas velocity and pressure gradient for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_idx.count,
        interp_ng, MFInfo(), *(m_eb->factory()[lev]));

      int components_count(0);

      // Copy fluid velocity
      MultiFab::Copy(*interp_ptr, *a_vel[lev], 0, interp_idx.vel, 3, interp_ng);
      components_count += interp_idx.ncomp_vel;

      // Copy pressure gradient
      MultiFab::Copy(*interp_ptr, gp_tmp, 0, interp_idx.grad_p, 3, interp_ng);
      components_count += interp_idx.ncomp_grad_p;

      // Copy fluid temperature
      if (include_Tf) {
        MultiFab::Copy(*interp_ptr, *a_T[lev], 0, interp_idx.Tf, 1, interp_ng);
        components_count += interp_idx.ncomp_Tf;
      }

      // Copy divtau
      if (include_divtau) {
        MultiFab::Copy(*interp_ptr, *a_divtau[lev], 0, interp_idx.divtau, 3, interp_ng);
        components_count += interp_idx.ncomp_divtau;
      }

      // Copy DufDt
      if (include_vm) {
        MultiFab::Copy(*interp_ptr, *a_DufDt[lev], 0, interp_idx.DufDt, 3, interp_ng);
        components_count += interp_idx.ncomp_DufDt;
      }

      AMREX_ALWAYS_ASSERT(interp_idx.count == components_count);

    } else {

      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_idx.count,
        interp_ng, MFInfo(), *(m_eb->particle_factory()[lev]));

      int components_count(0);

      // Copy fluid velocity
      interp_ptr->ParallelCopy(*a_vel[lev], 0, interp_idx.vel, 3, interp_ng, interp_ng);
      components_count += interp_idx.ncomp_vel;

      // Copy pressure gradient
      interp_ptr->ParallelCopy(gp_tmp, 0, interp_idx.grad_p, 3, interp_ng, interp_ng);
      components_count += interp_idx.ncomp_grad_p;

      // Copy fluid temperature
      if (include_Tf) {
        interp_ptr->ParallelCopy(*a_T[lev], 0, interp_idx.Tf, 1, interp_ng, interp_ng);
        components_count += interp_idx.ncomp_Tf;
      }

      // Copy divtau
      if (include_divtau) {
        interp_ptr->ParallelCopy(*a_divtau[lev], 0, interp_idx.divtau, 3, interp_ng, interp_ng);
        components_count += interp_idx.ncomp_divtau;
      }

      // Copy DufDt
      if (include_vm) {
        interp_ptr->ParallelCopy(*a_DufDt[lev], 0, interp_idx.DufDt, 3, interp_ng, interp_ng);
        components_count += interp_idx.ncomp_DufDt;
      }

      AMREX_ALWAYS_ASSERT(interp_idx.count == components_count);
    }

    // FillBoundary on interpolation MultiFab
    interp_ptr->FillBoundary(geom[lev].periodicity());

    BL_PROFILE_VAR("particle_deposition", particle_deposition);

    {
      const auto dx  = geom[lev].CellSizeArray();
      const auto dxi = geom[lev].InvCellSizeArray();
      const auto plo = geom[lev].ProbLoArray();

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(interp_ptr->Factory());

      const auto cellcent = &(factory.getCentroid());
      const auto bndrycent = &(factory.getBndryCent());
      const auto areafrac = factory.getAreaFrac();

      auto& plev = pc->GetParticles(lev);

      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {

        MFIXParticleContainer::PairIndex index(pti.index(), pti.LocalTileIndex());

        auto& ptile = plev[index];
        auto ptile_data = ptile.getParticleTileData();

        auto& particles = pti.GetArrayOfStructs();
        MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();

        const int np = particles.size();

        Box bx = pti.tilebox ();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  interp_fab = static_cast<EBFArrayBox const&>((*interp_ptr)[pti]);
        const EBCellFlagFab&  flags = interp_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {

          const auto& interp_array = interp_ptr->array(pti);

          const auto& flags_array = flags.array();

          const int grown_bx_is_regular = (flags.getType(amrex::grow(bx,1)) == FabType::regular);

          const Array4<const Real> empty_array;

          // Cell centroids
          const auto& ccent_fab = grown_bx_is_regular ? empty_array : cellcent->const_array(pti);
          // Centroid of EB
          const auto& bcent_fab = grown_bx_is_regular ? empty_array : bndrycent->const_array(pti);
          // Area fractions
          const auto& apx_fab = grown_bx_is_regular ? empty_array : areafrac[0]->const_array(pti);
          const auto& apy_fab = grown_bx_is_regular ? empty_array : areafrac[1]->const_array(pti);
          const auto& apz_fab = grown_bx_is_regular ? empty_array : areafrac[2]->const_array(pti);

          // We need this until we remove static attribute from mfix::gp0;
          const RealVect gp0_dev(gp0);

          int const explicit_update( (m_dem.solve() && (m_dem.implicit_drag()==0)) ? 1 : 0);

          int const idx_vel( interp_idx.vel );
          int const idx_gp( interp_idx.grad_p );
          int const idx_Tf( interp_idx.Tf );
          int const idx_divtau( interp_idx.divtau );
          int const idx_DufDt( interp_idx.DufDt );
          int const idx_count( interp_idx.count );

          int const idx_pc_vm_coeff(pc->m_runtimeRealData.vm_coeff);
          int const idx_pc_dudt(pc->m_runtimeRealData.acceleration);

          int const idx_pc_conv_coeff( pc->m_runtimeRealData.conv_coeff);
          int const idx_pc_energy_src( pc->m_runtimeRealData.energy_source);

          bool grown_bx_has_porous_media = (m_porous_media.nregions() > 0) ?
            m_porous_media.intersects_box(lev,amrex::grow(bx,1)) : false;

          const auto& pm_regions = grown_bx_has_porous_media ?
            m_porous_media.regions<run_on>() : PorousMediaRegions();

          ParallelFor(np, [pstruct,p_realarray,ptile_data,interp_array,gp0_dev,
            plo,dx,dxi,flags_array,ccent_fab,bcent_fab,apx_fab,apy_fab,apz_fab,
            grown_bx_is_regular,include_Tf,include_divtau,include_vm,explicit_update,
            dt=a_dt,idx_count,idx_vel,idx_gp,idx_Tf,idx_divtau,idx_DufDt,
            idx_pc_vm_coeff,idx_pc_dudt, idx_pc_conv_coeff, idx_pc_energy_src,
            grown_bx_has_porous_media,pm_regions,lev]
            AMREX_GPU_DEVICE (int p_id) noexcept
          {
            MFIXParticleContainer::ParticleType& particle = pstruct[p_id];

            // Local array storing interpolated values
            GpuArray<Real, interp_array_size> interp_loc;
            interp_loc.fill(0.);

            if (grown_bx_is_regular && !(grown_bx_has_porous_media)) {

              trilinear_interp(particle.pos(), interp_loc.data(),
                               interp_array, plo, dxi, idx_count);

            } else { // FAB not all regular or has porous media

              // Cell containing particle centroid
              const int iloc = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
              const int jloc = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
              const int kloc = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

              // The particle is in a covered cell or in a porous media
              if (flags_array(iloc,jloc,kloc).isCovered() ||
                  (grown_bx_has_porous_media && pm_regions.contains(lev,iloc,jloc,kloc)))
              {
                p_realarray[SoArealData::vel_source_x][p_id] = 0.0;
                p_realarray[SoArealData::vel_source_y][p_id] = 0.0;
                p_realarray[SoArealData::vel_source_z][p_id] = 0.0;

                if (include_vm) { ptile_data.m_runtime_rdata[idx_pc_vm_coeff][p_id] = 0.0; }

                if (include_Tf) {
                  ptile_data.m_runtime_rdata[idx_pc_conv_coeff][p_id] = 0.0;
                  ptile_data.m_runtime_rdata[idx_pc_energy_src][p_id] = 0.0;
                }

                return;

              // Cut or regular cell and none of the cells in the stencil is covered
              // (Note we can't assume regular cell has no covered cells in the stencil
              //      because of the diagonal case)
              } else {

                // Upper cell in trilinear stencil
                const int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5));
                const int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5));
                const int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5));

                // All cells in the stencil are regular. Use
                // traditional trilinear interpolation
                if (flags_array(i-1,j-1,k-1).isRegular() &&
                    flags_array(i  ,j-1,k-1).isRegular() &&
                    flags_array(i-1,j  ,k-1).isRegular() &&
                    flags_array(i  ,j  ,k-1).isRegular() &&
                    flags_array(i-1,j-1,k  ).isRegular() &&
                    flags_array(i  ,j-1,k  ).isRegular() &&
                    flags_array(i-1,j  ,k  ).isRegular() &&
                    flags_array(i  ,j  ,k  ).isRegular()) {

                  trilinear_interp(particle.pos(), interp_loc.data(),
                                   interp_array, plo, dxi, idx_count);

                // At least one of the cells in the stencil is cut or covered
                } else {
#if 0
                  // TODO: This was initially split for variables that may have known
                  // EB values (e.g., no-slip velocity). However, the results changed
                  // more than expected so now EB values are not used.
                  {
                    const int srccomp = 0;
                    const int dstcomp = 0;
                    const int numcomp = 3;

                    shepard_interp_eb(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                      flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                      interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
                  }
                  {
                    const int srccomp = 3;
                    const int dstcomp = 3;
                    const int numcomp = idx_count-3;

                    shepard_interp(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                   flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                   interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
                  }
#else
                  const int srccomp = 0;
                  const int dstcomp = 0;
                  const int numcomp = idx_count;

                  shepard_interp(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                 flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                 interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
#endif
                } // Cut cell
              } // Not covered
            } // if box not all regular

            Real pbeta = p_realarray[SoArealData::drag_coeff][p_id];
            Real pvol = SoArealData::volume(p_realarray[SoArealData::radius][p_id]);

            // Particle drag calculation.  We multiply the particle velocity
            // by "explicit_update" so that DEM uses the slip velocity. For PIC we
            // only want the fluid velocity as it uses a pseudo implicit
            // slip velocity for parcels.

            RealVect slip_vel(interp_loc[idx_vel], interp_loc[idx_vel+1], interp_loc[idx_vel+2]);
            RealVect gp_g(interp_loc[idx_gp], interp_loc[idx_gp+1], interp_loc[idx_gp+2]);

            if (explicit_update) {
              slip_vel[0] -= p_realarray[SoArealData::velx][p_id];
              slip_vel[1] -= p_realarray[SoArealData::vely][p_id];
              slip_vel[2] -= p_realarray[SoArealData::velz][p_id];
            }

            p_realarray[SoArealData::vel_source_x][p_id] =
              pbeta * slip_vel[0] - (gp_g[0] + gp0_dev[0]) * pvol;

            p_realarray[SoArealData::vel_source_y][p_id] =
              pbeta * slip_vel[1] - (gp_g[1] + gp0_dev[1]) * pvol;

            p_realarray[SoArealData::vel_source_z][p_id] =
              pbeta * slip_vel[2] - (gp_g[2] + gp0_dev[2]) * pvol;

            if (include_divtau) {
              p_realarray[SoArealData::vel_source_x][p_id] += interp_loc[idx_divtau  ]*pvol;
              p_realarray[SoArealData::vel_source_y][p_id] += interp_loc[idx_divtau+1]*pvol;
              p_realarray[SoArealData::vel_source_z][p_id] += interp_loc[idx_divtau+2]*pvol;
            }

            if (include_vm) { // Virtual mass force

              Real rhoCVp = ptile_data.m_runtime_rdata[idx_pc_vm_coeff][p_id];

              RealVect slip_acc(interp_loc[idx_DufDt], interp_loc[idx_DufDt+1], interp_loc[idx_DufDt+2]);

              if (explicit_update) {
                slip_acc[0] -= ptile_data.m_runtime_rdata[idx_pc_dudt  ][p_id];
                slip_acc[1] -= ptile_data.m_runtime_rdata[idx_pc_dudt+1][p_id];
                slip_acc[2] -= ptile_data.m_runtime_rdata[idx_pc_dudt+2][p_id];
              }

              p_realarray[SoArealData::vel_source_x][p_id] += rhoCVp*slip_acc[0];
              p_realarray[SoArealData::vel_source_y][p_id] += rhoCVp*slip_acc[1];
              p_realarray[SoArealData::vel_source_z][p_id] += rhoCVp*slip_acc[2];
            }

            if (include_Tf) {

              // gamma == (heat transfer coeff) * (particle surface area)
              Real const gamma = ptile_data.m_runtime_rdata[idx_pc_conv_coeff][p_id];
              Real const Tf = interp_loc[idx_Tf];

              ptile_data.m_runtime_rdata[idx_pc_energy_src][p_id] =
                gamma*(Tf - p_realarray[SoArealData::temperature][p_id]);

              //p_realarray[SoArealData::convection][p_id] =
              //  pgamma * (interp_loc[idx_Tf] - p_realarray[SoArealData::temperature][p_id]);
            }

          }); // particle loop
        } // FAB not covered
      } // pti
    } // omp region

    BL_PROFILE_VAR_STOP(particle_deposition);

    delete interp_ptr;

  } // lev

  // Reset velocity Dirichlet bc's to face values
  for (int lev(0); lev < nlev(); lev++) {
    bcs().set_velocity_bcs(lev, a_time, a_vel[lev]);
  }
}
