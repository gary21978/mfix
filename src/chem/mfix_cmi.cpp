#include <mfix_cmi.H>
#include <mfix_fluid.H>
#include <mfix_pc.H>
#include <mfix_mf_helpers.H>
#include <mfix_indexes_aux.H>
#include <mfix_material_interface.H>
#include <mfix_run_on.H>
#include <mfix_usr_reactions_rates_K.H>

#include <mfix_integrator.H>
#include <fe_integrator.H>
#include <be_integrator.H>
#include <vode_integrator.H>


using namespace amrex;


// Perform the whole chemistry update
void ChemistryManagementInterface::
calc_chemistry (Vector< Geometry > const& geom,
                Vector< MultiFab* > const& chem_txfr_out,
                Vector< MultiFab* > const& interp_ptr,
                const Real& therm_p_in,
                const Real time,
                const Real dt,
                const int verbose)
{
  BL_PROFILE("ChemistryManagementInterface::calc_chemistry()");

#ifdef OPTIMIZED_USR_RATES
  if (!rates_file_already_checked) {
    if (rates_aux::test_fluid_indexes(*m_fluid)) {
      reporter::Log(reporter::Error,__FILE__, __LINE__) << " Error:"
        << "The fluid species in the inputs file do not match the ones in the user reaction rates file!\n"
        << "Please correct the input deck or update the user reaction rates file.\n"
        << "To update the user reaction rates file run the python script in tools/Chemistry";
    }
    if (rates_aux::test_solids_indexes(*m_solids)) {
      reporter::Log(reporter::Error,__FILE__, __LINE__) << " Error:"
        << "The solids species in the inputs file do not match the ones in the user reaction rates file!\n"
        << "Please correct the input deck or update the user reaction rates file.\n"
        << "To update the user reaction rates file run the python script in tools/Chemistry";
    }
    if (rates_aux::test_lagrangian_indexes(*m_reactions)) {
      reporter::Log(reporter::Error,__FILE__, __LINE__) << " Error:"
        << "The chemical reactions in the inputs file do not match the ones in the user reaction rates file!\n"
        << "Please correct the input deck or update the user reaction rates file.\n"
        << "To update the user reaction rates file run the python script in tools/Chemistry";
    }
    if (rates_aux::test_eulerian_indexes(*m_reactions)) {
      reporter::Log(reporter::Error,__FILE__, __LINE__) << " Error:"
        << "The chemical reactions in the inputs file do not match the ones in the user reaction rates file!\n"
        << "Please correct the input deck or update the user reaction rates file.\n"
        << "To update the user reaction rates file run the python script in tools/Chemistry";
    }

    rates_file_already_checked = 1;
  }
#endif

#ifdef AMREX_DEBUG
  // Flag to check whether we're solving for the solids phase or not
  const int solve_solids = int(m_dem->solve() || m_pic->solve());
  // Flag to check whether we're solving lagrangian reactions
  const int lagrangian_reactions = (m_reactions->lagrangian_nreactions() > 0)? 1 : 0;
  // Flag to check whether we're solving eulerian reactions
  const int eulerian_reactions = (m_reactions->eulerian_nreactions() > 0)? 1 : 0;

  // Sanity check: if we're here we must be solving either the eulerian or the
  // lagrangian reactions
  AMREX_ASSERT(eulerian_reactions || lagrangian_reactions);

  // Sanity check: if we're here we are exclusively either not be solving for
  // eulerian reactions or, if we are solving for eulerian reactions, we must
  // also be solving for fluid, fluid density, fluid enthalpy and fluid species
  AMREX_ASSERT((!eulerian_reactions) ^
      (eulerian_reactions && m_fluid->solve() && m_fluid->solve_density() &&
       m_fluid->solve_enthalpy() && m_fluid->solve_species()));

  // Sanity check: if we're here we are exclusively either not be solving for
  // lagrangian reactions or, if we are solving for lagrangian reactions, we
  // must also be solving for solids, solids enthalpy and solids species
  AMREX_ASSERT((!lagrangian_reactions) ^
      (lagrangian_reactions && solve_solids &&
       m_solids->solve_enthalpy() && m_solids->solve_species()));
#endif

  // Save the start time for timing purposes
  const Real strttime = ParallelDescriptor::second();
  // Get the number of levels
  const int nlev = interp_ptr.size();

  // Perform the chemistry update only if we're solving chemical reactions
  if (m_reactions->solve()) {

    // Reset to zero the fluid chemistry transfer quantities
    for (int lev(0); lev < nlev; ++lev) {
      chem_txfr_out[lev]->setVal(0.);
    }

    // MultiFabs to store the number of particles in each fluid cell
    Vector<iMultiFab*> parts_nb_in_cell(nlev, nullptr);
    // MultiFabs to map each fluid cell to the first address in the vector
    // containing the particles IDs
    Vector<iMultiFab*> parts_idxs_in_cell(nlev, nullptr);
    // Vector containing the particles IDs ordered and grouped by the fluid
    // cells
    Vector<std::map<MFIXParticleContainer::PairIndex,
                    Gpu::DeviceVector<int>>> parts_idxs_vectorized(nlev);

    // Auxiliary MultiFabs that coincide with chem_txfr_out in case of same
    // grids between fluid and solids, or gets allocated on the particles grids
    // in case of different grids between fluid and solids
    Vector<MultiFab*> chem_txfr_mf(nlev, nullptr);

    // Get a vector of flags (one for each level) to discriminate whether the
    // fluid phas eand the particles have the same grids
    Vector<int> OnSameGrids(nlev, 1);

    // Loop over the levels
    for (int lev(0); lev < nlev; ++lev) {

      const DistributionMapping& dm( interp_ptr[lev]->DistributionMap() );
      const BoxArray&            ba( interp_ptr[lev]->boxArray() );

      OnSameGrids[lev] = int ((ba == chem_txfr_out[lev]->boxArray())   &&
                              (dm == chem_txfr_out[lev]->DistributionMap()));

      // If fluid and solids have the same grids
      if (OnSameGrids[lev]) {

        // Simply copy chem_txfr_out pointers into the auxiliary variables
        chem_txfr_mf[lev] = chem_txfr_out[lev];

      } else { // if fluid and solids have different grids

        // Allocate the auxiliary MultiFabs on the particles grid
        chem_txfr_mf[lev] = new MultiFab(ba, dm, chem_txfr_out[lev]->nComp(), 0,
          MFInfo(), interp_ptr[lev]->Factory());

        // Initialize the auxiliary MultiFabs to zero
        chem_txfr_mf[lev]->setVal(0.);
      }

      // If solving for particles
      if (m_dem->solve() || m_pic->solve()) {

        parts_nb_in_cell[lev] = new iMultiFab(ba, dm, 1, 0);
        parts_nb_in_cell[lev]->setVal(0);

        parts_idxs_in_cell[lev] = new iMultiFab(ba, dm, 1, 0);
        parts_idxs_in_cell[lev]->setVal(0);

        // Compute the particles information for each fluid cell
        calc_particles_within_cells(time, geom, interp_ptr, parts_nb_in_cell,
            parts_idxs_in_cell, parts_idxs_vectorized);
      }

    }

    // Perform the chemistry update
    calc_chem_txfr(time, dt, geom, chem_txfr_mf, interp_ptr, therm_p_in,
        parts_nb_in_cell, parts_idxs_in_cell, parts_idxs_vectorized);

    // Loop over the levels
    for (int lev(0); lev < nlev; ++lev) {
      // If fluid and solids are not on the same grids
      if (!OnSameGrids[lev]) {
        // Copy the chemistry transfer quantities from the solids grids to the
        // fluid grids
        chem_txfr_out[lev]->ParallelCopy(*chem_txfr_mf[lev], 0, 0,
            chem_txfr_mf[lev]->nComp(), 0, 0);
      }
    }

    // Free memory
    for (int lev(0); lev < nlev; ++lev) {
      if (!OnSameGrids[lev]) {
        delete chem_txfr_mf[lev];
      }
    }

    // Free memory
    if (m_dem->solve() || m_pic->solve()) {
      for (int lev(0); lev < nlev; ++lev) {
        delete parts_nb_in_cell[lev];
        delete parts_idxs_in_cell[lev];
      }
    }
  }

  // If the verbosity is greater than 1 print out the elapsed time
  if (verbose > 1) {

    // Get current time to compute the elapsed time since this function's start
    Real stoptime = ParallelDescriptor::second() - strttime;

    // Do parallel max reduction over the elapsed time
    ParallelDescriptor::ReduceRealMax(stoptime, ParallelDescriptor::IOProcessorNumber());

    // Master thread prints the elapsed time
    amrex::Print() << "Calc chemistry time: " << stoptime << '\n';
  }
}


// Compute the chemistry update
void ChemistryManagementInterface::
calc_chem_txfr (const Real time,
                const Real dt,
                Vector< Geometry > const& geom,
                Vector< MultiFab* > const& txfr_out,
                Vector< MultiFab* > const& interp_ptr,
                Real const& therm_p_in,
                Vector< iMultiFab* > const& parts_nb_in_cells,
                Vector< iMultiFab* > const& parts_idxs_in_cells,
                Vector< std::map< MFIXParticleContainer::PairIndex,
                                  Gpu::DeviceVector<int> > >& parts_idxs_vectorized)
{
  // Get the integration information for switching between different integrators
  const auto& integrator_parms = m_integrator_inputs.m_params;
  const int integrator_type = m_integrator_inputs.m_type;

  // If the integrator type is MFIXIntegrator
  if (integrator_type == IntegratorType::ForwardEuler) {
    // Call the ODE integration with the MFIXIntegrator
    calc_chem_txfr(time, dt, geom, txfr_out, interp_ptr, therm_p_in,
        parts_nb_in_cells, parts_idxs_in_cells, parts_idxs_vectorized,
        MFIXIntegrator(integrator_parms));

    // else if the integrator type is the StiffSolver::ForwardEuler
  } else if (integrator_type == IntegratorType::StiffSolver::FE) {
    // Call the ODE integration with the StiffSolver::ForwardEuler
    calc_chem_txfr(time, dt, geom, txfr_out, interp_ptr, therm_p_in,
        parts_nb_in_cells, parts_idxs_in_cells, parts_idxs_vectorized,
        ForwardEuler(integrator_parms));
  } else if (integrator_type == IntegratorType::StiffSolver::BE) {
    // Call the ODE integration with the StiffSolver::BackwardEuler
    calc_chem_txfr(time, dt, geom, txfr_out, interp_ptr, therm_p_in,
        parts_nb_in_cells, parts_idxs_in_cells, parts_idxs_vectorized,
        BackwardEuler(integrator_parms));
  } else if (integrator_type == IntegratorType::StiffSolver::VODE) {
    // Call the ODE integration with the StiffSolver::VODE
    calc_chem_txfr(time, dt, geom, txfr_out, interp_ptr, therm_p_in,
        parts_nb_in_cells, parts_idxs_in_cells, parts_idxs_vectorized,
        VODE(integrator_parms));
  }
}


// Compute the chemistry update
template <class F1>
void ChemistryManagementInterface::
calc_chem_txfr (const Real time,
                const Real dt,
                Vector< Geometry > const& geom,
                Vector< MultiFab* > const& txfr_out,
                Vector< MultiFab* > const& interp_ptr,
                Real const& therm_p_in,
                Vector< iMultiFab* > const& parts_nb_in_cells,
                Vector< iMultiFab* > const& parts_idxs_in_cells,
                Vector< std::map< MFIXParticleContainer::PairIndex,
                                  Gpu::DeviceVector<int> > >& parts_idxs_vectorized,
                const F1& integrator)
{
  using PairIndex = MFIXParticleContainer::PairIndex;

  BL_PROFILE("ChemistryManagementInterface::calc_chem_txfr()");

  // Get the number of levels
  const int nlev = interp_ptr.size();

  // Get fluid cell edge size and compute fluid cell volume
  const auto dx = geom[0].CellSizeArray();
  const auto reg_cell_vol = dx[0]*dx[1]*dx[2];

  // Get the eulerian and lagrangian reactions data
  const int lagrangian_reactions = (m_reactions->lagrangian_nreactions() > 0)? 1 : 0;
  const int eulerian_reactions = (m_reactions->eulerian_nreactions() > 0)? 1 : 0;

  // Set the maximum number of reactions as the maximum value between the
  // eulerian and lagrangian reactions. This is a strategy to save memory as
  // we will use the same memory space for both types of reactions
  const int max_nreactions = amrex::max(m_reactions->lagrangian_nreactions(),
                                        m_reactions->eulerian_nreactions());

  // Flags to check whether there are any fluid or solids species involved in
  // the chemical reactions. Start by setting the flags' value depending on the
  // eulerian and lagrangian flags
  short reacting_fluid_species = eulerian_reactions? 1 : 0;
  short reacting_solids_species = lagrangian_reactions? 1 : 0;

  // If any fluid species is found in any lagrangian reaction reactants or
  // products, then set the corresponding flag to true
  for (int n(0); (n < m_reactions->lagrangian_nreactions()) && (!reacting_fluid_species); ++n) {
    auto const& lagrangian_phases = m_reactions->get_lagrangian(n)->get_phases();
    for (auto phase: lagrangian_phases) {
      if (phase == ChemicalPhase::Fluid) {
        reacting_fluid_species = 1;
        break;
      }
    }
  }

  // If any solids species is found in any lagrangian reaction reactants or
  // products, then set the corresponding flag to true
  for (int n(0); (n < m_reactions->lagrangian_nreactions()) && (!reacting_solids_species); ++n) {
    auto const& lagrangian_phases = m_reactions->get_lagrangian(n)->get_phases();
    for (auto phase: lagrangian_phases) {
      if (phase == ChemicalPhase::Solid) {
        reacting_solids_species = 1;
        break;
      }
    }
  }

  // Flag to set whether we're solving for the fluid phase in the reactions
  const int solve_fluid = int(m_fluid->solve() && reacting_fluid_species);
  // Flag to set whether we're solving for the solids phase in the reactions
  const int solve_solids = int((m_dem->solve() || m_pic->solve()) && reacting_solids_species);

  // Get the number of fluid species
  const int nspecies_g = solve_fluid? m_fluid->nspecies() : 0;

  // Get the number of solids species
  const int nspecies_s = solve_solids? m_solids->nspecies() : 0;

  // Get the indexes to access the proper component for each different transfer
  // quantity for the fluid phase due to chemical reactions
  InterphaseChemTxfrIndexes chem_txfr_idxs(nspecies_g, m_reactions->solve());

  // Get the indexes to access the proper component in the MultiFab that
  // contains the fluid MultiFabs copied to the particles grids in case of
  // different fluid/solids grids
  DualGridAuxIndexes dualgrid_idxs(m_fluid->solve_enthalpy(), m_fluid->nspecies(),
      m_fluid->isMixture());

  // Get the fluid host/device parameters
  const auto& fluid_parms = m_fluid->parameters<run_on>();
  const auto fluid_props = m_fluid->props.data<run_on>();

  // Get the solids host/device parameters
  MFIXSolidsParms dummy_solids_parms;
  const auto& solids_parms = solve_solids? m_solids->parameters<run_on>() : dummy_solids_parms;
  const auto solids_props = m_solids->props.data<run_on>();

  // Get the indexes to access the NeighborParticleContainer runtime variables
  const auto runtime_indexes = solve_solids ? m_pc->m_runtimeRealData : runtimeRealData();

  // Get the eulerian and lagrangian reactions host/device parameters
  const auto& eulerian_reactions_parms = m_reactions->eulerian_parameters<run_on>();
  const auto& lagrangian_reactions_parms = m_reactions->lagrangian_parameters<run_on>();

  const int solids_update_type = m_integrator_inputs.m_solids_update_type;
  auto solids_thresholds = m_integrator_inputs.m_solids_thresholds;

  //***************************************************************************
  //
  //***************************************************************************

  const int report_mass_balance = m_report_mass_balance;

  if (report_mass_balance) {
    m_fluid_mass_change.clear();
    m_solids_mass_change.clear();

    m_fluid_mass_change.resize(nspecies_g, 0.);
    m_solids_mass_change.resize(nspecies_s, 0.);
  }

  Real* fluid_mass_change_ptr = report_mass_balance? m_fluid_mass_change.dataPtr() : nullptr;
  Real* solids_mass_change_ptr = report_mass_balance? m_solids_mass_change.dataPtr() : nullptr;

  //***************************************************************************
  //
  //***************************************************************************

  // Loop over the levels
  for (int lev = 0; lev < nlev; lev++) {

    // Loop over the Boxes
    for (MFIter mfi(*txfr_out[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      // Get the particles' data
      PairIndex index(mfi.index(), mfi.LocalTileIndex());

      MFIXParticleContainer::ParticleTileType dummy_ptile;
      auto& ptile = solve_solids? m_pc->GetParticles(lev)[index] : dummy_ptile;

      //Access to added variables
      auto ptile_data = ptile.getParticleTileData();

      // Get the current Box
      Box bx = mfi.tilebox();

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(txfr_out[lev]->Factory());
      const auto& flags = factory.getMultiEBCellFlagFab();

      // If the current Box is not fully covered
      if (flags[mfi].getType(bx) != FabType::covered) {

        // Get the number of cells in this box
        const int bx_numPts = bx.numPts();

        // Get the Array4-s for the fluid chemistry transfer quantities, the
        // flags, and the volume fraction
        const auto& txfr_array = txfr_out[lev]->array(mfi);
        const auto& flags_array = flags.const_array(mfi);
        const auto& volfrac_array = (factory.getVolFrac()).const_array(mfi);

        const int fluid_is_mixture = m_fluid->isMixture();

        // Dummy variables needed in case we do not solve for solids or the
        // fluid is not a mixture of species
        const Array4<const Real> real_dummy_arr;

        // Get the fluid Array4-s for each fluid variable
        const auto& vel_g_array = solve_fluid ?
          interp_ptr[lev]->const_array(mfi, dualgrid_idxs.vel_g) : real_dummy_arr;

        const auto& ep_g_array = solve_fluid ?
          interp_ptr[lev]->const_array(mfi, dualgrid_idxs.ep_g) : real_dummy_arr;

        const auto& ro_g_array = solve_fluid ?
          interp_ptr[lev]->const_array(mfi, dualgrid_idxs.ro_g) : real_dummy_arr;

        const auto& T_g_array = solve_fluid ?
          interp_ptr[lev]->const_array(mfi, dualgrid_idxs.T_g) : real_dummy_arr;

        const auto& h_g_array = solve_fluid ?
          interp_ptr[lev]->const_array(mfi, dualgrid_idxs.h_g) : real_dummy_arr;

        const auto& X_gk_array = (solve_fluid && fluid_is_mixture) ?
          interp_ptr[lev]->const_array(mfi, dualgrid_idxs.X_gk) : real_dummy_arr;

        // Get the fluid thermodynamic pressure at this level
        const Real thermo_p = therm_p_in;

        // Dummy variables needed in case we do not solve for solids or the
        // fluid is not a mixture of species
        amrex::Array4<int> int_dummy_arr;

        // Get the Array4-s from the MultiFabs that contain the information
        // about which particles are in which fluid cells
        const auto& parts_nb_array = solve_solids ?
          parts_nb_in_cells[lev]->const_array(mfi) : int_dummy_arr;

        const auto& parts_idxs_array = solve_solids ?
          parts_idxs_in_cells[lev]->const_array(mfi) : int_dummy_arr;

        int* parts_idxs_ptr(nullptr);
        if (solve_solids) {
          Gpu::DeviceVector<int>& parts_idxs_vect = parts_idxs_vectorized[lev][index];
          parts_idxs_ptr = parts_idxs_vect.dataPtr();
        }

        // Define the memory handler class, which allocates the host/device
        // memory for the chemistry ODE integration
        ChemMemoryHandler mem_handler(bx_numPts, max_nreactions, nspecies_g,
            nspecies_s, solve_fluid, solve_solids, lagrangian_reactions, integrator);

        // Get the pointers to the allocated host/device memory
        Real* rates_mem = mem_handler.rates_mem();
        Real* variables_mem = mem_handler.variables_mem();
        Real* rhs_mem = mem_handler.rhs_mem();
        Real* y_old_mem = mem_handler.y_old_mem();

        // Get the flag if we're using device shared memory or not
        const int use_shared_mem = mem_handler.use_shared_mem();

        // Get the number of variables for the ODE integrator auxiliary
        // variables for the three different types: short, 1D Real, 2D Real
        const int naux_short1D = mem_handler.naux_short1D();
        const int naux_real1D = mem_handler.naux_real1D();
        const int naux_real2D = mem_handler.naux_real2D();

        // Get the pointers to the allocated ODE integrator auxiliary variables
        // host/device memory
        short* aux_short1D_mem = (naux_short1D > 0)? mem_handler.aux_short1D_mem() : nullptr;
        Real* aux_real1D_mem = (naux_real1D > 0)? mem_handler.aux_real1D_mem() : nullptr;
        Real* aux_real2D_mem = (naux_real2D > 0)? mem_handler.aux_real2D_mem() : nullptr;

        // Depending if we're doing GPU or non-GPU compilation, we need to
        // launch the ParallelFor loop in different ways because of the
        // possibility to use GPU shared memory which is defined differently for
        // SYCL, and cannot be called for non-GPU compilations.
#ifdef AMREX_USE_GPU
        // In case of GPU compilation, the amrex::launch routine allows to
        // specify the GPU shared memory size
        amrex::launch(mem_handler.nblocks(), mem_handler.nthreads_per_block(),
            mem_handler.needed_sm_bytes(), Gpu::gpuStream(),
#else
        // In case of non-GPU compilation we can simply call the
        // amrex::ParallelFor routine
        amrex::ParallelFor(bx_numPts,
#endif
            [ptile_data,parts_idxs_ptr,parts_nb_array,parts_idxs_array,
             flags_array,ep_g_array,ro_g_array,T_g_array,vel_g_array,X_gk_array,
             nspecies_g,thermo_p,nspecies_s,lagrangian_reactions,eulerian_reactions,
             txfr_array,fluid_parms,fluid_props,solids_parms,solids_props,reg_cell_vol,bx,bx_numPts,
             max_nreactions,solve_solids,volfrac_array,dt,runtime_indexes,
             chem_txfr_idxs,eulerian_reactions_parms,time,solve_fluid,
             lagrangian_reactions_parms,h_g_array,integrator,rates_mem,y_old_mem,
             variables_mem,rhs_mem,use_shared_mem,aux_short1D_mem,aux_real1D_mem,
             aux_real2D_mem,naux_short1D,naux_real1D,naux_real2D,fluid_mass_change_ptr,
             solids_mass_change_ptr,report_mass_balance,solids_update_type,
             solids_thresholds]
#ifdef AMREX_USE_GPU
          // non-SYCL, GPU compilation
          AMREX_GPU_DEVICE () noexcept
        {
          // Allocate and define the material interface class for the current
          // cell
          MaterialInterface mic;
          mic.setup(bx_numPts, max_nreactions, nspecies_g, nspecies_s, solve_fluid, solve_solids,
              lagrangian_reactions, use_shared_mem, naux_short1D, naux_real1D, naux_real2D);
          Gpu::SharedMemory<Real> gsm; // the GPU shared memory
          mic.setup_memory(gsm, rates_mem, variables_mem, rhs_mem, y_old_mem, aux_short1D_mem,
              aux_real1D_mem, aux_real2D_mem);
#else
          // non-GPU compilation
          (int idx) noexcept
        {
          // Allocate and define the material interface class for the current
          // cell
          MaterialInterface mic;
          mic.setup(bx_numPts, max_nreactions, nspecies_g, nspecies_s, solve_fluid, solve_solids,
              lagrangian_reactions, use_shared_mem, naux_short1D, naux_real1D, naux_real2D);
          mic.setup_memory(idx, rates_mem, variables_mem, rhs_mem, y_old_mem, aux_short1D_mem,
              aux_real1D_mem, aux_real2D_mem);
#endif

          // Get the global index from the memory handler. For GPU runs, this
          // index is the current GPU thread global index; for non-GPU runs,
          // this index is equal to the ParallelFor loop index
          const int glob_idx = mic.memory_handler().glob_idx();

          // We check that the global index is within the limits to avoid
          // accesses out of memory. This is needed in GPU runs as the number of
          // threads launched can be larger than the size of the problem
          if (glob_idx < bx_numPts) {

            // Retrieve the 3D indexes ijk from the 1D index glob_idx (remember
            // that we have vectorized the problem over the current Box)
            IntVect ijk = get_local_indexes(glob_idx, bx);

            // If the current fluid cell is not covered
            if (!flags_array(ijk).isCovered()) {

              // Allocate and define the ODE integrator for the current fluid
              // cell
              F1 gpu_integrator = integrator;
              gpu_integrator.define(mic.memory_handler(), ijk, time);
              gpu_integrator.setup_aux_memory();

              // Get the fluid cell volume and volume fraction
              const Real cell_volume = reg_cell_vol;
              const Real vfrac = volfrac_array(ijk);

              // Allocate and define the eulerian reactions data
              MFIXReactionsData eulerian_reactions_data(mic.memory_handler(),
                  eulerian_reactions_parms, time, dt, cell_volume, vfrac);
              // Define the material interface class pointer to the eulerian
              // reactions class
              mic.setup_eulerian_reactions(&eulerian_reactions_data,
                  eulerian_reactions);

              // Allocate and define the lagrangian reactions data
              MFIXReactionsData lagrangian_reactions_data(mic.memory_handler(),
                  lagrangian_reactions_parms, time, dt, cell_volume, vfrac);
              // Define the material interface class pointer to the lagrangian
              // reactions class
              mic.setup_lagrangian_reactions(&lagrangian_reactions_data,
                  lagrangian_reactions);

              // Allocate and define the fluid data
              MFIXFluidData fluid_data(mic.memory_handler(), fluid_parms, fluid_props, ijk,
                  ep_g_array, ro_g_array, vel_g_array, T_g_array, h_g_array,
                  X_gk_array, thermo_p, cell_volume, vfrac, txfr_array,
                  chem_txfr_idxs, report_mass_balance, fluid_mass_change_ptr);
              // Define the material interface class pointer to the fluid data
              // class
              mic.setup_fluid(&fluid_data, solve_fluid);

              // Allocate and define the solids data
              MFIXSolidsData solids_data(mic.memory_handler(), solids_parms, solids_props, ijk,
                  parts_nb_array, parts_idxs_array, parts_idxs_ptr, ptile_data,
                  runtime_indexes, solids_update_type, solids_thresholds,
                  gpu_integrator.integration_status(), report_mass_balance,
                  solids_mass_change_ptr);
              // Define the material interface class pointer to the solids data
              // class
              mic.setup_solids(&solids_data, solve_solids);

              // Compute the chemical reactions update by integrating the
              // chemistry ODE system
              mic.compute_reactions(&gpu_integrator, dt);

            } // end if cell is not covered
          } // end if glob_idx < numPts

        }); // end ParallelFor

        // If we have allocated device global memory, synchronize host and
        // device to avoid the deallocation of some temporary data and memory
        // that has been allocated in this scope for the chemistry computation
        if (mem_handler.glob_mem_size() != 0)
          Gpu::synchronize();

      } // end if entire FAB not covered

      if (solve_solids) {
        removeInvalidParticles(ptile);
      }

    } // end mfi loop

  } // end loop over the nlev levels
}
