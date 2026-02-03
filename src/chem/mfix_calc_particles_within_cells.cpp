#include <AMReX_EBFArrayBox.H>

#include <mfix_cmi.H>
#include <mfix_mf_helpers.H>
#include <mfix_indexes_aux.H>


using namespace amrex;


// Compute the information regarding which particles are contained in which
// fluid cells.
void ChemistryManagementInterface::
calc_particles_within_cells (const Real /*time*/,
                             Vector< Geometry > const& geom,
                             Vector< MultiFab* > const& interp_ptr,
                             Vector< iMultiFab* > const& mf_particles_nb_per_cell,
                             Vector< iMultiFab* > const& mf_global_indxs_per_cell,
                             Vector< std::map< MFIXParticleContainer::PairIndex,
                                               Gpu::DeviceVector<int> > >& global_indxs_vectorized)
{
  // Define some useful aliases
  using PairIndex = MFIXParticleContainer::PairIndex;
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::calc_particles_within_cells()");

  // Get the number of levels
  const int nlev = interp_ptr.size();

  // Loop over the levels
  for (int lev(0); lev < nlev; ++lev) {

    // Get a reference to this level vector of particles' IDs
    std::map<PairIndex,
             Gpu::DeviceVector<int>>& indxs_vectorized = global_indxs_vectorized[lev];

    // Loop over the boxes
    for (MFIXParIter pti(*m_pc, lev); pti.isValid(); ++pti)
    {
      // Insert a new element key:value in the map for this level vector of
      // particles' IDs
      PairIndex index(pti.index(), pti.LocalTileIndex());
      indxs_vectorized[index] = Gpu::DeviceVector<int>();
    }

    // Get the inverse of the fluid cell edge size
    const auto dxi = geom[lev].InvCellSizeArray();
    // Get the lower limits of the geometry
    const auto plo = geom[lev].ProbLoArray();

    // Loop over the boxes
    for (MFIXParIter pti(*m_pc, lev); pti.isValid(); ++pti) {

      // Get the particles' data
      PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& particles = pti.GetArrayOfStructs();
      MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();
      const int np = particles.size();

      // Get current Box
      const Box& bx = pti.tilebox();

      // This is to check efficiently if this tile contains any eb stuff
      const EBFArrayBox& fab = static_cast<EBFArrayBox const&>((*interp_ptr[lev])[pti]);
      const EBCellFlagFab& flags = fab.getEBCellFlagFab();

      // Check if this Box is fully covered
      if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {

        // Local variables for storing the particles-per-cell information. Note
        // that we use vectors (instead of Array4-s) because later we'll call
        // the amrex::inclusive_scan function on these data structures.
        Gpu::DeviceVector<Long> particles_nb_per_cell;
        Gpu::DeviceVector<Long> global_indxs_per_cell;

#ifdef AMREX_DEBUG
        // Local variables to double-check the correctness of the computation in
        // case of DEBUG mode compilation
        Gpu::DeviceScalar<Long> particles_in_covered_cells(0);
        Long* particles_in_covered_cells_ptr = particles_in_covered_cells.dataPtr();
#endif

        // Get the number of cells in this Box
        const int bx_numPts = bx.numPts();

        // Initialize the local variables
        particles_nb_per_cell.clear();
        particles_nb_per_cell.resize(bx_numPts, 0);

        global_indxs_per_cell.clear();
        global_indxs_per_cell.resize(bx_numPts, 0);

        // Get a pointer to the vector of the number of particles for each fluid
        // cell
        Long* particles_nb_per_cell_ptr = particles_nb_per_cell.dataPtr();

        // Get the Array4 for the EB flags
        const auto& flags_array = flags.const_array();

        // Loop over the particles
        amrex::ParallelFor(np, [pstruct,plo,dxi,bx,particles_nb_per_cell_ptr,
#ifdef AMREX_DEBUG
            flags_array,particles_in_covered_cells_ptr]
#else
            flags_array]
#endif
          AMREX_GPU_DEVICE (int p_id) noexcept
        {
          // Get a reference to the current particle
          MFIXParticleContainer::ParticleType& particle = pstruct[p_id];

          // Compute the fluid cell indices where the particle is located
          const int i = int(Math::floor((particle.pos(0) - plo[0])*dxi[0]));
          const int j = int(Math::floor((particle.pos(1) - plo[1])*dxi[1]));
          const int k = int(Math::floor((particle.pos(2) - plo[2])*dxi[2]));

          // If the fluid cell is not covered
          if (!flags_array(i,j,k).isCovered()) {

            // Convert the 3D indexes ijk into a 1D index
            const Long n = get_global_index(IntVect(i,j,k), bx);

            // Add 1 atomically for this particle to the fluid cell where the
            // particle is located
            HostDevice::Atomic::Add(&particles_nb_per_cell_ptr[n], Long(1));

#ifdef AMREX_DEBUG
          } else {
            // For debug purposes, in case a particle is in  a covered cell, add
            // 1 to the covered cell location
            HostDevice::Atomic::Add(particles_in_covered_cells_ptr, Long(1));
#endif
          }

        }); // end ParallelFor loop

        // Compute the global indices for each fluid cell using the inclusive
        // scan function from AMReX. The resulting vector of global indices for
        // each fluid cell is important as we will use it as a map from the
        // fluid cell to the vector containing the particles IDs ordered and
        // grouped by fluid cells.
        Gpu::inclusive_scan(particles_nb_per_cell.begin(),
            particles_nb_per_cell.end(), global_indxs_per_cell.begin());

#ifdef AMREX_DEBUG
        // Sanity check: the total number of particles in the fluid cells
        // computed above must match the total number of particles from the
        // MFIXParticleContainer object
        {
          Long local_np(0);
          const int vec_size = bx_numPts;

#ifdef AMREX_USE_GPU
          // If the data structures are on the GPU, we need to copy that data
          // back to the host to perform the sanity check
          Gpu::HostVector<Long> global_indxs_per_cell_host(vec_size, 0);
          Gpu::copy(Gpu::deviceToHost, global_indxs_per_cell.begin(),
              global_indxs_per_cell.end(), global_indxs_per_cell_host.begin());

          local_np = global_indxs_per_cell_host[vec_size-1];
#else
          local_np = global_indxs_per_cell[vec_size-1];
#endif

          Long covered_np = particles_in_covered_cells.dataValue();
          // Sanity check
          AMREX_ALWAYS_ASSERT((local_np+covered_np) == np);
        }
#endif

        // Here we copy the quantities computed above from the local
        // DeviceVectors to the MultiFabs that will be used later on in the
        // chemistry computation
        Long* global_indxs_per_cell_ptr = global_indxs_per_cell.dataPtr();

        const Array4<int>& particles_per_cell_array = mf_particles_nb_per_cell[lev]->array(pti);
        const Array4<int>& global_indxs_array = mf_global_indxs_per_cell[lev]->array(pti);

        // Loop over the current Box
        amrex::ParallelFor(bx, [plo,dxi,global_indxs_array,flags_array,bx,
            particles_per_cell_array,global_indxs_per_cell_ptr,
            particles_nb_per_cell_ptr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          // Check whether this cell is covered
          if (!flags_array(i,j,k).isCovered()) {

            // Get the 1D index from the 3D indices for the current cell of
            // this Box
            const Long n = get_global_index(IntVect(i,j,k), bx);

            // Copy from the vectors to the MultiFabs
            particles_per_cell_array(i,j,k) = particles_nb_per_cell_ptr[n];
            global_indxs_array(i,j,k) = global_indxs_per_cell_ptr[n];

          } else {
            // Do nothing.
          }

        }); // end ParallelFor loop

        // Reset the data structure where we store the particles IDs for each
        // fluid cell
        indxs_vectorized[index].clear();
        indxs_vectorized[index].resize(np, -1);
        int* indxs_vectorized_ptr = indxs_vectorized[index].dataPtr();

        // Loop over the particles
        amrex::ParallelFor(np, [pstruct,plo,dxi,global_indxs_array,flags_array,
            particles_nb_per_cell_ptr,global_indxs_per_cell_ptr,bx,indxs_vectorized_ptr]
          AMREX_GPU_DEVICE (int p_id) noexcept
        {
        // Get a reference to this particle
          MFIXParticleContainer::ParticleType& particle = pstruct[p_id];

          // Compute the indices of the cell where current particle is located
          const int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
          const int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
          const int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

          // Check if this cell is not covered
          if (!flags_array(i,j,k).isCovered()) {

            // Get the 1D index from the 3D indices of this cell in the Box
            const Long n = get_global_index(IntVect(i,j,k), bx);
            // Get the global index of the particles contained in this cell.
            // This global index is computed through the inclusive scan
            // operation performed above, and represents how many particles
            // there are in the non-covered fluid cells up to this cell included
            const Long idx = global_indxs_per_cell_ptr[n];
            // Atomically subtract 1 from the number of particles in this cell.
            // Also return the number of particles in this cell before the
            // subtraction has happened. In this way we can compute the global
            // index for this particle in the order of the fluid cells.
            const int nb_temp = Gpu::Atomic::Add(&particles_nb_per_cell_ptr[n], Long(-1));

            // Compute the global index for this particle in the order of the
            // fluid cells.
            const Long nn = idx - Long(nb_temp);
            // Assign this particle's ID to the position associated to this cell
            indxs_vectorized_ptr[nn] = p_id;

          } else {
            // Do nothing. In this case, we might have to throw an error
          }

        }); // end ParallelFor loop
      } // end if entire FAB not covered

      // host/device synchronization
      Gpu::synchronize();

    } // en loop over the boxes
  } // end the loop over the levels lev

}
