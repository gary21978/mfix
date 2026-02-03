#include <mfix_pc.H>
#include <mfix_des_K.H>
#include <mfix_pc_interactions_K.H>
#include <mfix_des_roll_friction_K.H>
#include <mfix_pc_updates_K.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_calc_cell.H>

using namespace amrex;
using namespace Solvers;

void MFIXParticleContainer::
EvolveParticles (int const lev, Real const dt,
                 RealVect& gravity,
                 EBFArrayBoxFactory* ebfactory,
                 EBFArrayBoxFactory* particle_ebfactory,
                 const MultiFab* ls_phi,
                 const int ls_refinement,
                 LoadBalance *const a_loadbalance,
                 const MultiFab* a_T_eb,
                 int& nsubsteps)
{
  if (m_dem.rolling_friction() == RollingFrictionModel::ModelA) {

    EvolveParticles(lev, dt, gravity, ebfactory, particle_ebfactory,
                    ls_phi, ls_refinement, a_loadbalance, a_T_eb,
                    nsubsteps, ModelA(m_dem.rolling_friction_coeff()));

  } else if (m_dem.rolling_friction() == RollingFrictionModel::ModelB) {

    EvolveParticles(lev, dt, gravity, ebfactory, particle_ebfactory,
                    ls_phi, ls_refinement, a_loadbalance, a_T_eb,
                    nsubsteps, ModelB(m_dem.rolling_friction_coeff()));

  } else if (m_dem.rolling_friction() == RollingFrictionModel::None) {

    EvolveParticles(lev, dt, gravity, ebfactory, particle_ebfactory,
                    ls_phi, ls_refinement, a_loadbalance, a_T_eb,
                    nsubsteps, NoRollingFriction());

  } else {

    Abort("Unknown Rolling Friction model");
  }
}

template <class F1>
void MFIXParticleContainer::
EvolveParticles ( int const lev, Real const dt,
                  RealVect& gravity,
                  EBFArrayBoxFactory* ebfactory,
                  EBFArrayBoxFactory* particle_ebfactory,
                  const MultiFab* ls_phi,
                  const int ls_refinement,
                  LoadBalance* const a_loadbalance,
                  const MultiFab* a_T_eb,
                  int& nsubsteps,
                  F1 RollingFriction)
{
    BL_PROFILE_REGION_START("mfix_dem::EvolveParticles()");
    BL_PROFILE("mfix_dem::EvolveParticles()");

    amrex::Print() << "Evolving particles on level: " << lev
                   << " ... with fluid dt " << dt << std::endl;


    /****************************************************************************
     * Geometry                                                                 *
     ***************************************************************************/

    const Real* dx = Geom(lev).CellSize();

    const int has_T_eb = bcs().eb_parms().has_temperature();
    if (has_T_eb) { AMREX_ASSERT( a_T_eb != nullptr ); }

    /****************************************************************************
     * Init substeps                                                            *
     ***************************************************************************/

    Real subdt;
    // des_init_time_loop(&dt, &nsubsteps, &subdt);
    if ( dt >= m_dem.dtsolid() )
    {
       nsubsteps = static_cast<int>(amrex::Math::ceil(dt / static_cast<amrex::Real>(m_dem.dtsolid())));
       subdt     = dt / nsubsteps;
    } else {
       nsubsteps = 1;
       subdt     = dt;
    }

    /****************************************************************************
     * Get particle EB geometric info
     ***************************************************************************/
    const FabArray<EBCellFlagFab>* flags = &(particle_ebfactory->getMultiEBCellFlagFab());


    /****************************************************************************
     * Init temporary storage:                                                  *
     *   -> particle-particle, and particle-wall forces                         *
     *   -> particle-particle, and particle-wall torques                        *
     ***************************************************************************/
    std::map<PairIndex, Gpu::DeviceVector<Real>> tow;
    std::map<PairIndex, Gpu::DeviceVector<Real>> fc;

    std::map<PairIndex, Gpu::DeviceVector<Real>> cond;

    std::map<PairIndex, bool> tile_has_walls;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const Box& bx = pti.tilebox();

        PairIndex index(pti.index(), pti.LocalTileIndex());

        tow[index]  = Gpu::DeviceVector<Real>();
        fc[index]   = Gpu::DeviceVector<Real>();
        cond[index] = Gpu::DeviceVector<Real>();

        // Determine if this particle tile actually has any walls
        bool has_wall = false;

        if ((particle_ebfactory != NULL)
            && ((*flags)[pti].getType(amrex::grow(bx,1)) == FabType::singlevalued))
        {
            has_wall = true;
        }
        else
        {
            // We need this test for the case of an inflow boundary:
            // inflow does not appear in the EBFactory but
            // the particles see it as a wall

            // Create the nodal refined box based on the current particle tile
            Box refined_box(amrex::convert(amrex::refine(bx,ls_refinement), IntVect{1,1,1}));

            // Set tol to 1/2 dx
            Real tol = amrex::min(dx[0], amrex::min(dx[1], dx[2])) / 2;

            Real ls_min_over_box = ((*ls_phi)[pti]).min<run_on>(refined_box,0);

            if (ls_min_over_box < tol) has_wall = true;
        }

        tile_has_walls[index] = has_wall;
    }

    /****************************************************************************
     * Iterate over sub-steps                                                   *
     ***************************************************************************/

    // Particle inflow
    if (ebfactory != NULL) {
      mfix_pc_inflow(lev, 1, 0, dt, m_solids.solve_enthalpy(), ebfactory);
    }

    const Real abstol = newton_abstol;
    const Real reltol = newton_reltol;
    const int maxiter = newton_maxiter;

    int n = 0; // Counts sub-steps
    while (n < nsubsteps)
    {
        // Redistribute particles ever so often BUT always update the neighbour
        // list (Note that this fills the neighbour list after every
        // redistribute operation)
        if (m_solids.update_momentum()) {
          if (n % 25 == 0) {
              clearNeighbors();
              Redistribute(0, 0, 0, 0, false); // Do not remove negatives
              fillNeighbors();
#ifdef AMREX_USE_GPU
              if (reduceGhostParticles) {
                selectActualNeighbors(MFIXCheckPair(m_dem.neighborhood()));
                updateNeighbors(true);
              }
#endif
              if( m_dem.pneig_flag() ) {
#if MFIX_POLYDISPERSE
                  buildNeighborList(MFIXCheckPolyPair(SoAintData::ptype,m_dem.nptypes(),m_dem.pneighdata()),
                                    SoAintData::ptype, m_dem.prefratdata(), m_dem.nptypes(), false);
#else
                  amrex::Abort("MFIX not built with POLYDISPERSE support");
#endif
              } else {
                  buildNeighborList(MFIXCheckPair(m_dem.neighborhood()), false);
              }
          } else {
              updateNeighbors();
          }
        }

        /********************************************************************
         * Particles routines                                               *
         *******************************************************************/
        BL_PROFILE_VAR("particles_computation", particles_computation);
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // Timer used for load-balancing
            Real const timer_start( ParallelDescriptor::second() );

            //const Box& bx = pti.tilebox(); // UNUSED_VARIABLE
            PairIndex index(pti.index(), pti.LocalTileIndex());

            auto& plev  = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            const auto ntp = aos.size();
            const int  nrp = GetParticles(lev)[index].numRealParticles();

            // For multi-grid neighbor list search, we must
            // loop over the ghost particles to find all coll pairs.
            const int  nlp = (m_dem.pneig_flag()) ? ntp : nrp;

            auto& soa = pti.GetStructOfArrays();
            auto p_realarray = soa.realarray();
            auto p_intarray  = soa.intarray();

            //Access to added variables
            auto ptile_data = ptile.getParticleTileData();

            // Number of particles including neighbor particles
            int ntot = nrp;

            // Particle-particle (and particle-wall) forces and torques. We need
            // these to be zero every time we start a new batch (i.e tile and
            // substep) of particles.
            tow[index].clear();
            tow[index].resize(3*ntot, 0.0);

            Real* tow_ptr  = tow[index].dataPtr();

            fc[index].clear();
            fc[index].resize(3*ntot, 0.0);

            Real* fc_ptr   = fc[index].dataPtr();

            cond[index].clear();
            cond[index].resize(ntot, 0.0);

            Real* cond_ptr = cond[index].dataPtr();

            auto& geom = this->Geom(lev);
            GpuArray<int,3> periodic = geom.isPeriodicArray();

            if (m_solids.update_momentum()) {

              // For debugging: keep track of particle-particle (pfor) and
              // particle-wall (wfor) forces

              // BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);

              const auto dxi = geom.InvCellSizeArray();
              const auto plo = geom.ProbLoArray();
              const auto& phiarr = ls_phi->array(pti);

              const auto& T_eb = (has_T_eb) ? a_T_eb->const_array(pti)
                                            : Array4<Real const>{};

              const int walls_in_tile = tile_has_walls[index];

              auto nbor_data = m_neighbor_list[lev][index].data();

              const auto& solids_parms = m_solids.parameters<run_on>();
              const int solve_enthalpy = m_solids.solve_enthalpy();

              const int pft_neighbor_idx = m_runtimeRealData.pft_neighbor_idx;
              const int max_contacts_tan_history = m_runtimeRealData.max_contacts;
              const int cpu_id_idx = m_runtimeIntData.cpu_id_idx;
              const int pid_idx = m_runtimeIntData.pid_idx;
              const int wall_pid_idx = m_runtimeIntData.wall_pid_idx;

              Gpu::DeviceVector<bool> touch(nlp*max_contacts_tan_history, false);
              bool* touch_ptr = touch.dataPtr();

              int const include_conduction(solve_enthalpy && solids_parms.get_do_pfp_cond());

              // now we loop over the neighbor list and compute the forces
              amrex::ParallelFor(nlp,
                  [nrp,pstruct,p_realarray,p_intarray,fc_ptr,tow_ptr,cond_ptr,
                   nbor_data,subdt,ntot,walls_in_tile,ls_refinement,phiarr,plo,
                   dxi,solids_parms,solve_enthalpy, has_T_eb, T_eb, include_conduction,
                   RollingFriction,mew=m_dem.mew(),mew_w=m_dem.mew_w(),kn=m_dem.kn(),
                   kn_w=m_dem.kn_w(),etan=m_dem.etan(),etan_w=m_dem.etan_w(),
                   k_g=m_dem.k_g_dem(),kt=m_dem.kt(),etat=m_dem.etat(),
                   kt_w=m_dem.kt_w(),etat_w=m_dem.etat_w(),
                   tan_history=m_dem.tan_history(),ptile_data,max_contacts_tan_history,
                   pft_neighbor_idx,pid_idx,cpu_id_idx,wall_pid_idx,touch_ptr]
                AMREX_GPU_DEVICE (int i) noexcept
                {
                    auto particle = pstruct[i];

                    RealVect pos1(particle.pos());

                    RealVect total_force(0.);
                    RealVect total_tow_force(0.);

                    const int istate(p_intarray[SoAintData::state][i]);

                    //**********************************************************
                    // Particle-wall collisions
                    //**********************************************************
                    if (walls_in_tile && (i < nrp)) {

                      int const include_wall_cond( include_conduction && has_T_eb );

                      particle_walls(p_realarray, p_intarray, i,
                          ls_refinement, phiarr, plo, dxi, subdt, pos1,
                          solids_parms, include_wall_cond, T_eb, k_g,
                          kn_w, kt_w, etan_w, etat_w, mew_w, total_force, total_tow_force,
                          cond_ptr, istate, RollingFriction,
                          tan_history, ptile_data, max_contacts_tan_history,
                          pft_neighbor_idx, cpu_id_idx, pid_idx, wall_pid_idx,
                          touch_ptr);

                    } // tile has walls

                    //**********************************************************
                    // Particle-particle collisions
                    //**********************************************************
                    const auto neighbs = nbor_data.getNeighbors(i);

                    particle_particles(particle, neighbs, p_realarray,
                        p_intarray, i, solve_enthalpy, subdt, nrp,
                        solids_parms, k_g, kn, kt, etan, etat, mew, total_force,
                        total_tow_force, fc_ptr, tow_ptr, cond_ptr,
                        ntot, istate, RollingFriction,
                        tan_history, ptile_data,
                        max_contacts_tan_history, pft_neighbor_idx, cpu_id_idx,
                        pid_idx, wall_pid_idx, touch_ptr);
              });

              Gpu::Device::synchronize();

              // BL_PROFILE_VAR_STOP(calc_particle_collisions);

              // BL_PROFILE_VAR("des::solve_particle_velocity_and_position()", des_time_march);

            }

            /********************************************************************
             * Move particles based on collision forces and torques             *
             *******************************************************************/

            const auto p_lo = Geom(lev).ProbLoArray();
            const auto p_hi = Geom(lev).ProbHiArray();

            GpuArray<int,6> lo_hi_bc{m_boundary_conditions.domain_bc(0),
                                     m_boundary_conditions.domain_bc(1),
                                     m_boundary_conditions.domain_bc(2),
                                     m_boundary_conditions.domain_bc(3),
                                     m_boundary_conditions.domain_bc(4),
                                     m_boundary_conditions.domain_bc(5)};

            int const include_vm( m_runtimeRealData.contains_vm() );
            int const idx_pc_vm_coeff( m_runtimeRealData.vm_coeff );

            int const include_acc( m_runtimeRealData.contains_acc() );
            int const idx_pc_acc( m_runtimeRealData.acceleration );

            const int nspecies_s = m_solids.nspecies();

            const int idx_X_sn = m_runtimeRealData.X_sn;

            const int idx_conv_coeff = m_runtimeRealData.conv_coeff;
            const int idx_energy_src = m_runtimeRealData.energy_source;

            const int update_momentum = m_solids.update_momentum();
            const int solve_enthalpy = m_solids.solve_enthalpy();

            const Real enthalpy_source = m_solids.enthalpy_source();

            const int solid_is_a_mixture = m_solids.isMixture();
            const auto& solids_parms = m_solids.parameters<run_on>();
            const auto solids_props = m_solids.props.data<run_on>();

            bool implicit_drag = m_dem.implicit_drag();

            amrex::ParallelFor(nrp, [pstruct,p_realarray,p_intarray,subdt,
                ptile_data,nspecies_s,idx_X_sn,fc_ptr,cond_ptr,p_hi,p_lo,ntot,
                lo_hi_bc,periodic,enthalpy_source,update_momentum,gravity,tow_ptr,
                solid_is_a_mixture,solids_parms,solids_props,solve_enthalpy,implicit_drag,
                abstol,reltol,maxiter, include_vm, idx_pc_vm_coeff, include_acc,
                idx_pc_acc, idx_conv_coeff, idx_energy_src, dt ]
              AMREX_GPU_DEVICE (int i) noexcept
            {
              ParticleType& p = pstruct[i];

              //***************************************************************
              // Second step: update particles' positions and velocities
              //***************************************************************
              if (update_momentum) {

                part_momentum_update(p, ptile_data, p_realarray, p_intarray, i, subdt, ntot,
                    fc_ptr, tow_ptr, gravity, p_lo, p_hi, lo_hi_bc, periodic, implicit_drag, dt,
                    include_vm, idx_pc_vm_coeff, include_acc, idx_pc_acc);
              }

              //***************************************************************
              // Third step: update particles' temperature
              //***************************************************************
              if (solve_enthalpy) {

                part_enthalpy_update(ptile_data, p_realarray, i, idx_X_sn,
                    idx_conv_coeff, idx_energy_src, solids_props, subdt, cond_ptr,
                    enthalpy_source, abstol, reltol, maxiter, 1);

              }
            });

            Gpu::synchronize();

            usr2_des(nrp, ptile);

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
        }
        BL_PROFILE_VAR_STOP(particles_computation);

        // Update substep count
        n += 1;

    } // end of loop over substeps

    amrex::Print() << "done. \n";

    BL_PROFILE_REGION_STOP("mfix_dem::EvolveParticles()");
}


void MFIXParticleContainer::
PostEvolveParticles (int const lev)
{
    // Redistribute particles at the end of all substeps (note that the particle
    // neighbour list needs to be reset when redistributing).
    if (m_solids.update_momentum()) {

      clearNeighbors();
      Redistribute(0, 0, 0, 1, false); // Remove negatives

      auto& plev = GetParticles(lev);
      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto& ptile = plev[index];
        removeInvalidParticles(ptile);
      }
    }

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

        const int nrp   = pti.numParticles();
        void* particles = pti.GetArrayOfStructs().data();

        usr3_des(nrp,particles);
    }
}
