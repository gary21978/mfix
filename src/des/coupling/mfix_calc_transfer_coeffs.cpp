#include <mfix_coupling.H>

#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_shepard_K.H>
#include <mfix_mf_helpers.H>
#include <mfix_run_on.H>

#include <mfix_des_drag_K.H>
#include <mfix_des_added_mass_coeff_K.H>
#include <mfix_des_conv_coeff_K.H>

using namespace amrex;

void CouplingOp::
calc_transfer_coeffs ( MFIXFluidPhase const& a_fluid,
                       MFIXParticleContainer* a_pc,
                       GeometryData const& a_geomData )
{
  const auto& fluid_parms = a_fluid.parameters<run_on>();
  const auto fluid_props = a_fluid.props.data<run_on>();

  DragModel const* drag_model = getDragModel();

  if ( drag_model->WenYu() ) {

    calc_transfer_coeffs(a_fluid, a_pc, a_geomData,
      ComputeDragWenYu(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
        fluid_props, a_geomData));

  } else if ( drag_model->Gidaspow() ) {

    calc_transfer_coeffs(a_fluid, a_pc, a_geomData,
      ComputeDragGidaspow(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
        fluid_props, a_geomData));

  } else if ( drag_model->BVK2() ) {

    calc_transfer_coeffs(a_fluid, a_pc, a_geomData,
      ComputeDragBVK2(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
        fluid_props, a_geomData));

  } else if ( drag_model->SyamOBrien() ) {

    amrex::Real const c1( drag_model->SyamOBrien_coeff_c1() );
    amrex::Real const d1( drag_model->SyamOBrien_coeff_d1() );

    calc_transfer_coeffs(a_fluid, a_pc, a_geomData,
      ComputeDragSyamOBrien1988(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
        fluid_props, a_geomData, c1, d1));

  } else if ( drag_model->UserDrag() ) {

    calc_transfer_coeffs(a_fluid, a_pc, a_geomData,
      ComputeDragUser(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
        fluid_props, a_geomData));

  } else {

    amrex::Abort("Invalid drag type");

  }

}


template <typename F1>
void CouplingOp::
calc_transfer_coeffs ( MFIXFluidPhase const& a_fluid,
                       MFIXParticleContainer* a_pc,
                       GeometryData const& a_geomData,
                       F1 DragCoeff)
{
  const auto& fluid_parms = a_fluid.parameters<run_on>();
  const auto fluid_props = a_fluid.props.data<run_on>();

  if (include_virtual_mass()) {

    VirtualMassModel const* virtual_mass = getVirtualMassModel();

    if (virtual_mass->Null()) {

      // This is a 'do nothing' option that is primarily useful for testing.
      calc_transfer_coeffs(a_fluid, a_pc, a_geomData, DragCoeff,
        AddedMassCoeff_Null(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
          fluid_props, a_geomData));

    } else if (virtual_mass->Constant()) {

      calc_transfer_coeffs(a_fluid, a_pc, a_geomData, DragCoeff,
        AddedMassCoeff_Constant(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
          fluid_props, a_geomData, virtual_mass->getConstantCoeff()));

    } else if (virtual_mass->Zuber()) {

      calc_transfer_coeffs(a_fluid, a_pc, a_geomData, DragCoeff,
        AddedMassCoeff_Zuber(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
          fluid_props, a_geomData));

    } else if (virtual_mass->Nijssen()) {

      calc_transfer_coeffs(a_fluid, a_pc, a_geomData, DragCoeff,
        AddedMassCoeff_Nijssen(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
          fluid_props, a_geomData));

    } else {

      amrex::Abort("Invalid Virtual Mass model.");

    }
  } else {

    calc_transfer_coeffs(a_fluid, a_pc, a_geomData, DragCoeff,
      AddedMassCoeff_Null(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
        fluid_props, a_geomData));

  }
}



template <typename F1, typename F2>
void CouplingOp::
calc_transfer_coeffs ( MFIXFluidPhase const& a_fluid,
                       MFIXParticleContainer* a_pc,
                       GeometryData const& a_geomData,
                       F1 DragCoeff, F2 VirtualMassCoeff)
{
  const auto& fluid_parms = a_fluid.parameters<run_on>();
  const auto fluid_props = a_fluid.props.data<run_on>();

  if (include_energy()) {

    ConvectionModel const* convection_model = getConvectionModel();

    if ( convection_model->RanzMarshall() ) {
      calc_transfer_coeffs(a_fluid, a_pc, a_geomData, DragCoeff, VirtualMassCoeff,
        ComputeConvRanzMarshall(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
          fluid_props, a_geomData));

    } else if (convection_model->Gunn() ) {

      calc_transfer_coeffs(a_fluid, a_pc, a_geomData, DragCoeff, VirtualMassCoeff,
        ComputeConvGunn(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
          fluid_props, a_geomData));

    } else if ( convection_model->NullConvection() ) {

      calc_transfer_coeffs(a_fluid, a_pc, a_geomData, DragCoeff, VirtualMassCoeff,
        NullConvectionCoeff(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
          fluid_props, a_geomData));

    } else {

      amrex::Abort("Invalid convection type");

    }

  } else {

    calc_transfer_coeffs(a_fluid, a_pc, a_geomData, DragCoeff, VirtualMassCoeff,
      NullConvectionCoeff(a_pc->m_runtimeRealData, *m_interp_idxs, fluid_parms,
        fluid_props, a_geomData));

  }
}


template <typename F1, typename F2, typename F3>
void CouplingOp::
calc_transfer_coeffs ( MFIXFluidPhase const& a_fluid,
                       MFIXParticleContainer* a_pc,
                       GeometryData const& /*a_geomData*/,
                       F1 DragCoeff, F2 VirtualMassCoeff,
                       F3 ConvectionCoeff)
{
  using PairIndex = MFIXParticleContainer::PairIndex;
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("CouplingOp::calc_transfer_coeff()");

  // Set negative index to skip if virtual mass is not needed
  int const include_vm(a_pc->m_runtimeRealData.contains_vm());
  int const idx_pc_vm_coeff( a_pc->m_runtimeRealData.vm_coeff);

  if (include_vm) {
    AMREX_ASSERT( idx_pc_vm_coeff >= 0);
    AMREX_ASSERT( idx_pc_vm_coeff < a_pc->m_runtimeRealData.count );
  }

  // convection coefficient
  int const idx_pc_conv_coeff( a_pc->m_runtimeRealData.conv_coeff );

  const auto& fluid_parms = a_fluid.parameters<run_on>();

  int const idx_interp_vel(m_interp_idxs->vel_g);
  int const idx_interp_epf(m_interp_idxs->ep_g);
  int const idx_interp_rho(m_interp_idxs->ro_g);
  int const idx_interp_Tf(m_interp_idxs->T_g);
  int const idx_interp_Xfk(m_interp_idxs->X_gk);

  int const contains_Tf(include_energy());
  int const contains_Xf(include_species());
  int const nspecies(m_interp_idxs->nspecies);

  for (int lev = 0; lev < a_pc->numLevels(); lev++) {

    const auto dxi = a_pc->Geom(lev).InvCellSizeArray();
    const auto dx  = a_pc->Geom(lev).CellSizeArray();
    const auto plo = a_pc->Geom(lev).ProbLoArray();

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(m_interp[lev]->Factory());

    const auto cellcent = &(factory.getCentroid());
    const auto bndrycent = &(factory.getBndryCent());
    const auto areafrac = factory.getAreaFrac();

    auto& plev = a_pc->GetParticles(lev);

    for (MFIXParIter pti(*a_pc, lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& ptile = plev[index];
      auto ptile_data = ptile.getParticleTileData();

      auto& particles = pti.GetArrayOfStructs();
      MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();
      auto p_intarray = soa.intarray();

      const int np = particles.size();

      Box bx = pti.tilebox();

      // This is to check efficiently if this tile contains any eb stuff
      const EBFArrayBox& interp_fab = static_cast<EBFArrayBox const&>((*m_interp[lev])[pti]);
      const EBCellFlagFab& flags = interp_fab.getEBCellFlagFab();

      if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {

        const Array4<const Real> empty_array;

        const int interp_comp = m_interp[lev]->nComp();
        const auto& interp_array = m_interp[lev]->const_array(pti);

        const auto& flags_array  = flags.array();

        const int grown_bx_is_regular = (flags.getType(amrex::grow(bx,1)) == FabType::regular);

        // Cell centroids
        const auto& ccent_fab = grown_bx_is_regular ? empty_array : cellcent->const_array(pti);
        // Centroid of EB
        const auto& bcent_fab = grown_bx_is_regular ? empty_array : bndrycent->const_array(pti);
        // Area fractions
        const auto& apx_fab = grown_bx_is_regular ? empty_array : areafrac[0]->const_array(pti);
        const auto& apy_fab = grown_bx_is_regular ? empty_array : areafrac[1]->const_array(pti);
        const auto& apz_fab = grown_bx_is_regular ? empty_array : areafrac[2]->const_array(pti);

#ifdef AMREX_USE_GPU
        // Get from AMReX the maximum number of GPU threads in each block of
        // threads
        constexpr int nthreads_per_block = AMREX_GPU_MAX_THREADS;
        // Compute the number of blocks of threads needed to perform the whole
        // computation. Divide the number of particles by the number of threads
        // per block, and get the ceil of the result.
        int nblocks = static_cast<int>(amrex::Math::ceil(Real(np) / nthreads_per_block));

        // Compute the amount of shared memory (in bytes) needed by each block
        // of threads within the GPU kernel. Multiply the number of
        // interpolation components by the number of threads per block, by the
        // amount of bytes occupied by an object of type amrex::Real
        std::size_t needed_sm_bytes = interp_comp*nthreads_per_block*sizeof(Real);

        // Get the available amount (in bytes) of shared memory for the current
        // GPU device
        std::size_t available_sm_per_block = Gpu::Device::sharedMemPerBlock();

        // By default, try to use the GPU shared memory
        int use_shared_mem(1);

        // The GPU global memory
        Gpu::DeviceVector<Real> glob_mem;
        Real* glob_mem_ptr(nullptr);

        // If the needed shared memory is larger than the shared memory space
        // available on the current GPU device, use the GPU global memory
        if (needed_sm_bytes > available_sm_per_block) {
          // Set the shared memory flag to false
          use_shared_mem = 0;
          // Reset to zero the amount of needed shared memory
          needed_sm_bytes = 0;

          // Allocate the global memory. The size of needed memory is equal to
          // the number of interpolation components by the number of particles
          glob_mem.resize(interp_comp*np);
          // Assign the pointer to the GPU global memory
          glob_mem_ptr = glob_mem.dataPtr();
        }

        // Launch the GPU kernel
        amrex::launch(nblocks, nthreads_per_block, needed_sm_bytes, Gpu::gpuStream(),

#else

        // When not running on GPUs, shared memory flag is set to false and the
        // global memory pointer is set to nullptr
        int use_shared_mem(0);
        Real* glob_mem_ptr(nullptr);

        // Launch the loop over the particles
        amrex::ParallelFor(np,

#endif

          [pstruct,p_realarray,p_intarray,ptile_data,interp_array,
           DragCoeff, VirtualMassCoeff, ConvectionCoeff,
           plo,dxi,np,nspecies,interp_comp, fluid_parms,flags_array,grown_bx_is_regular,dx,apx_fab,
           apy_fab,apz_fab,ccent_fab,bcent_fab,use_shared_mem,glob_mem_ptr,
           idx_interp_vel, idx_interp_epf, idx_interp_rho, idx_interp_Tf, idx_interp_Xfk,
           contains_Tf, contains_Xf, include_vm, idx_pc_vm_coeff, idx_pc_conv_coeff]
#ifdef AMREX_USE_GPU

        AMREX_GPU_DEVICE () noexcept
        {
          // Get the thread global index
          const int p_id = blockDim.x*blockIdx.x+threadIdx.x;
          // Get the thread local index (local to its block of threads)
          const int thread_id = threadIdx.x;
          // Get the pointer to the GPU shared memory
          Gpu::SharedMemory<Real> shared_memory;
          Real* sm = shared_memory.dataPtr();


          Real* interp_local(nullptr);

          // Do work only if this thread global index is valid. This check is
          // needed because the number of GPU blocks of threads was computed as
          // the ceil of the division of the number of particles by the number
          // of threads in a block. Thus the global index of a thread might
          // exceed np. Threads with global index larger than np should not
          // perform any work, otherwise they might access invalid memory
          // addresses.
          if (p_id < np) {

            if (use_shared_mem) {

              // If shared memory is used, assign unique shared memory areas to
              // each thread. This is achieved by assigning the shared memory
              // space starting at address thread_id*interp_comp. Thus, each
              // thread will read/write in the shared memory space between
              // element thread_id*interp_comp and (thread_id+1)*interp_comp-1,
              // extremes included
              interp_local = &(sm[thread_id*interp_comp]);

            } else { // if not use_shared_mem

              // If global memory is used, assign unique global memory areas to
              // each thread. This is achieved by assigning the global memory
              // space starting at address p_id*interp_comp. Thus, each thread
              // will read/write in the global memory space between element
              // p_id*interp_comp and (p_id+1)*interp_comp-1, extremes included
              interp_local = &(glob_mem_ptr[p_id*interp_comp]);

            } // end if use_shared_mem

#else

        (int p_id) noexcept
        {
          // When not running on GPUs, allocate local memory and use it for
          // interpolation purposes
          Vector<Real> interp(interp_comp, 0.);
          Real* interp_local = interp.data();

          {

#endif

            MFIXParticleContainer::ParticleType& particle = pstruct[p_id];

            // Indices of cell where particle is located
            const int iloc = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
            const int jloc = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
            const int kloc = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

            if (grown_bx_is_regular) {

              trilinear_interp(particle.pos(), interp_local, interp_array,
                plo, dxi, interp_comp);

            } else { // FAB not all regular

              // No drag force for particles in covered cells.
              if (flags_array(iloc,jloc,kloc).isCovered()) {

                // drag variable
                p_realarray[SoArealData::drag_coeff][p_id] = 0.;

                // convection-related enthalpy txfr variable
                if (contains_Tf) {
                  ptile_data.m_runtime_rdata[idx_pc_conv_coeff ][p_id] = 0.0;
                }

                // virtual-mass coefficient
                if ( include_vm ) {
                  ptile_data.m_runtime_rdata[idx_pc_vm_coeff][p_id] = 0.0;
                }

                // Nothing else to do. Return
                return;

              // Cut or regular cell and none of the cells in the stencil is
              // covered (Note we can't assume regular cell has no covered
              // cells in the stencil because of the diagonal case)
              } else {

                // Upper cell in trilinear stencil
                int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5));
                int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5));
                int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5));

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

                  trilinear_interp(particle.pos(), interp_local, interp_array,
                      plo, dxi, interp_comp);

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
                    const int numcomp = interp_comp-3; // eps_f, rho_f, T_f, X_gk

                    shepard_interp(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                   flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                   interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
                  }
#else
                  const int srccomp = 0;
                  const int dstcomp = 0;
                  const int numcomp = interp_comp; // vel_g, eps_f, rho_f, T_f, p_g

                  shepard_interp(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                      flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                      interp_array, interp_local, srccomp, dstcomp, numcomp);

#endif
                } // Cut cell
              } // Not covered
            } // type of FAB


            Real vol_p = SoArealData::volume(p_realarray[SoArealData::radius][p_id]);

            //  interp_array, interp_local
            Real const drag_coeff = DragCoeff(p_id, particle, p_intarray, p_realarray,
                                              ptile_data, interp_array, interp_local);

            p_realarray[SoArealData::drag_coeff][p_id] = drag_coeff * vol_p;

            // virtual-mass coefficient
            if ( include_vm ) {

              Real const rho_f( interp_local[idx_interp_rho] );

              Real const vm_coeff = VirtualMassCoeff(p_id, particle, p_intarray, p_realarray,
                                                     ptile_data, interp_array, interp_local);

              ptile_data.m_runtime_rdata[idx_pc_vm_coeff][p_id] = rho_f * vm_coeff * vol_p;
            }

            // convection coefficient
            if(contains_Tf) {

              Real const rad_p( p_realarray[SoArealData::radius][p_id] );
              Real const Sa( 4.0*M_PI*rad_p*rad_p );

              Real const conv_coeff = ConvectionCoeff( p_id, particle, p_intarray, p_realarray,
                                                       ptile_data, interp_array, interp_local);

              ptile_data.m_runtime_rdata[idx_pc_conv_coeff ][p_id] = Sa * conv_coeff;

            } // solve enthalp

          }
        }); // pid
      } // if entire FAB not covered
    } // pti
  } // lev
}
