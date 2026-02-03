#include <mfix_reporter.H>
#include <mfix_run_on.H>
#include <mfix_calc_cell.H>

#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_ic.H>

#include <mfix_pc_generator.H>

using namespace amrex;

void MFIXParticleContainer::
InitParticlesAscii (const std::string& file)
{

  int const lev = 0;

  const RealVect dx(Geom(lev).CellSize());
  Real multi_particle_volume(0.);

  { Real const cell_vol = dx[0]*dx[1]*dx[2];

    Real const np_per_cell_at_pack = (m_pic.solve() ? m_pic.parcels_per_cell_at_pack() :
      ( m_dem.cg_dem() ? m_dem.cg_particles_per_cell_at_pack() : 1.));

    Real const ep_cp = (m_pic.solve() ? m_pic.ep_cp() : ( m_dem.cg_dem() ? 0.65 : 0.));

    multi_particle_volume = (cell_vol * ep_cp) / np_per_cell_at_pack;

  }

  // only read the file on the IO proc
  if (ParallelDescriptor::IOProcessor())
  {
    std::ifstream ifs;
    ifs.open(file, std::ios::in);

    if (!ifs.good())
      amrex::FileOpenFailed(file);

    int np = -1;
    ifs >> np >> std::ws;

    // Issue an error if nparticles = 0 is specified
    if ( np == -1 ){
      Abort("\nCannot read number of particles from particle_input.dat: file is corrupt.\
                   \nPerhaps you forgot to specify the number of particles on the first line?");
    }

    // we add all the particles to grid 0 and tile 0 and let
    // Redistribute() put them in the right places.
    const int grid = 0;
    const int tile = 0;

    auto& particles = DefineAndReturnParticleTile(lev,grid,tile);
    particles.resize(np);

    Gpu::HostVector<ParticleType> host_particles(np);

    std::array<Gpu::HostVector<Real>, SoArealData::count> host_realarrays;
    std::array<Gpu::HostVector<int>, SoAintData::count> host_intarrays;

    for (int comp(0); comp < SoArealData::count; ++comp)
      host_realarrays[comp].resize(np);

    for (int comp(0); comp < SoAintData::count; ++comp)
      host_intarrays[comp].resize(np);

    int  pstate, pphase;
    Real velx, vely, velz;
    Real pradius, pdensity, pomega;

    pstate = state::normal;

    int max_particle_phase(-1);

    for (int i = 0; i < np; i++)
    {
      // Read from input file
      ifs >> pphase;

      max_particle_phase = amrex::max(max_particle_phase, pphase);

      ifs >> host_particles[i].pos(0);
      ifs >> host_particles[i].pos(1);
      ifs >> host_particles[i].pos(2);
      ifs >> pradius;
      ifs >> pdensity;
      ifs >> velx;
      ifs >> vely;
      ifs >> velz;

      host_realarrays[SoArealData::velx][i]   = velx;
      host_realarrays[SoArealData::vely][i]   = vely;
      host_realarrays[SoArealData::velz][i]   = velz;

      // Compute other particle properties
      pomega = 0.0;

      // Set id and cpu for this particle
      host_particles[i].id()  = ParticleType::NextID();
      host_particles[i].cpu() = ParallelDescriptor::MyProc();

      // Set other particle properties
      host_intarrays[SoAintData::phase][i]        = pphase;
      host_intarrays[SoAintData::state][i]        = pstate;
#if MFIX_POLYDISPERSE
      host_intarrays[SoAintData::ptype][i]        = 0;
#endif

      host_realarrays[SoArealData::density][i]    = pdensity;
      host_realarrays[SoArealData::radius][i]     = pradius;
      host_realarrays[SoArealData::omegax][i]     = pomega;
      host_realarrays[SoArealData::omegay][i]     = pomega;
      host_realarrays[SoArealData::omegaz][i]     = pomega;

      // Compute the weight using the "multi particle volume" and
      // the "single particle" diameter.

      host_realarrays[SoArealData::statwt][i] = 1.;

      if (multi_particle_volume > 0.) {
        Real const dp = 2.*pradius;
        host_realarrays[SoArealData::statwt][i] =
          multi_particle_volume / ( (M_PI/6.0)*dp*dp*dp);
      }

      // Initialize these for I/O purposes
      host_realarrays[SoArealData::drag_coeff][i]    = 0.0;
      host_realarrays[SoArealData::vel_source_x][i]  = 0.0;
      host_realarrays[SoArealData::vel_source_y][i]  = 0.0;
      host_realarrays[SoArealData::vel_source_z][i]  = 0.0;
      host_realarrays[SoArealData::temperature][i]   = 0.0;

      if (!ifs.good())
          amrex::Abort("Error initializing particles from Ascii file. \n");
    }

    // NOTE : No need to do a ParallelDescriptor::ReduceIntMax on
    // max_particle_phase because we're reading the particle_input.dat only on
    // the IO proc

    if (m_dem.solve() && max_particle_phase > m_dem.NPHASE())
      amrex::Abort("One or more particle in the particle_input.dat has a phase number that is not present in the inputs file");
    else if (m_pic.solve() && max_particle_phase > m_pic.NPHASE())
      amrex::Abort("One or more particle in the particle_input.dat has a phase number that is not present in the inputs file");

    auto& aos = particles.GetArrayOfStructs();
    Gpu::DeviceVector<ParticleType>& gpu_particles = aos();

    // Copy particles from host to device
    Gpu::copyAsync(Gpu::hostToDevice, host_particles.begin(), host_particles.end(), gpu_particles.begin());

    auto& soa = particles.GetStructOfArrays();
    auto p_realarray = soa.realarray();
    auto p_intarray = soa.intarray();

    // Copy particles from host to device
    for (int comp(0); comp < SoArealData::count; ++comp) {
      Gpu::copyAsync(Gpu::hostToDevice, host_realarrays[comp].begin(),
          host_realarrays[comp].end(), &(p_realarray[comp][0]));
    }

    // Copy particles from host to device
    for (int comp(0); comp < SoAintData::count; ++comp) {
      Gpu::copyAsync(Gpu::hostToDevice, host_intarrays[comp].begin(),
          host_intarrays[comp].end(), &(p_intarray[comp][0]));
    }

    // Add components for each of the runtime variables
    const int start = SoArealData::count;
    for (int comp(0); comp < m_runtimeRealData.count; ++comp)
      particles.push_back_real(start+comp, np, 0.);
  }

  Redistribute(0,0,0,0,true);
}


void MFIXParticleContainer::
InitParticlesAuto (EBFArrayBoxFactory* particle_ebfactory)
{
  int lev = 0;

  const RealVect dx(Geom(lev).CellSize());
  const RealVect plo(Geom(lev).ProbLo());

  auto& ics = m_initial_conditions.ic();

  // Store particle count totals by IC region
  std::vector<long> total_np(ics.size(), 0);

  const Real tolerance = std::numeric_limits<Real>::epsilon();

  const auto& flags = particle_ebfactory->getMultiEBCellFlagFab();

  const int allow_ic_regions_overlap = m_initial_conditions.allow_overlap();

  Real multi_particle_volume(0.);

  { Real const cell_vol = dx[0]*dx[1]*dx[2];

    Real const np_per_cell_at_pack = (m_pic.solve() ? m_pic.parcels_per_cell_at_pack() :
      ( m_dem.cg_dem() ? m_dem.cg_particles_per_cell_at_pack() : 1.));

    Real const ep_cp = (m_pic.solve() ? m_pic.ep_cp() : ( m_dem.cg_dem() ? 0.65 : 0.));

    multi_particle_volume = (cell_vol * ep_cp) / np_per_cell_at_pack;
  }

  // double check if this goes out of MFIter loop
  MFIXICRegions ic_regions;

  for (int icv(0); icv < ics.size(); icv++) {

    auto& ic = ics[icv];
    auto& ic_solids = ic.solids;
    auto& ic_fluid = ic.fluid;

    AMREX_ALWAYS_ASSERT(ic_solids.size() <= 1);

    if (Math::abs(ic_fluid.volfrac-1) > tolerance) {

      const RealBox* ic_region = ic.region;
      const int ic_region_added = ic_regions.add(*ic_region, m_initial_conditions.allow_overlap());

      if (ic_region_added) {

        // Copy IC regions data from host to device
        const int ic_regions_size = ic_regions.size();
        Gpu::DeviceVector<RealBox> ic_regions_gpu(ic_regions_size);
        Gpu::copy(Gpu::hostToDevice, ic_regions.begin(), ic_regions.end(), ic_regions_gpu.begin());
        const RealBox* ic_regions_ptr = ic_regions_gpu.dataPtr();

        for (int lcs(0); lcs < ic_solids.size(); lcs++) {

          if (ic_solids[lcs].volfrac > tolerance) {

            const int phase = ic_solids[lcs].phase;

            ParticlesGenerator particles_generator( plo, dx,
                ic_regions_ptr, ic_regions_size, allow_ic_regions_overlap,
                *(m_initial_conditions.ic(icv).region),
                *(m_initial_conditions.ic(icv).get_solid(phase)),
                m_initial_conditions.ic(icv).packing,
                m_initial_conditions.has_granular_temperature(icv),
                m_dem.solve(), m_pic.solve(), m_dem.cg_dem(),
                multi_particle_volume);

            const Box ic_box = calc_ic_box(Geom(lev), ic_region);

            // This uses the particle tile size. Note that the default is to tile so if we
            //      remove the true and don't explicitly add false it will still tile
            for (MFIter mfi = MakeMFIter(lev,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

              const Box& tilebx = mfi.tilebox();
              auto& particles = DefineAndReturnParticleTile(lev, mfi);

              // Now that we know pcount, go ahead and create a particle container for this
              // grid and add the particles to it
              if(flags[mfi].getType(tilebx) != FabType::covered && tilebx.intersects(ic_box)) {

                const Box bx = tilebx & ic_box;

                const int id = ParticleType::NextID();
                const int cpu = ParallelDescriptor::MyProc();

                // This is particles in this grid for this IC region
                int pcount = particles_generator.generate(bx, id, cpu, particles);

                // Update the particles NextID
                ParticleType::NextID(id+pcount);

                // Add components for each of the runtime variables
                const int start = SoArealData::count;
                for (int comp(0); comp < m_runtimeRealData.count; ++comp)
                  particles.push_back_real(start+comp, pcount, 0.);

                total_np[icv] += static_cast<long>(pcount);
              } // if Fab is not covered and tilebox intersects IC region

              removeInvalidParticles(particles);

            } // MFIter loop

            break; // only one solid phase per icv is allowed
          } // ep_s > 0
        } // loop over solids
      } // if IC region was added
    } // ep_g < 1
  } // loop over ICs

  ParallelDescriptor::ReduceLongSum(total_np.data(), ics.size());

  long total_numparticle = 0;
  for (int icv(0); icv < ics.size(); icv++) {
    m_initial_conditions.set_particle_count(icv, total_np[icv]);
    total_numparticle += total_np[icv];
  }

  amrex::Print() << "Total number of generated particles: "
    << reporter::FormatWithCommas(total_numparticle) << std::endl;

  // We shouldn't need this if the particles are tiled with one tile per grid, but otherwise
  // we do need this to move particles from tile 0 to the correct tile.
  Redistribute(0,0,0,0,true);

  // We've already assigned a normal velocity with zero mean and a standard
  // deviation of 1 to all the particles we generated.
  if (m_initial_conditions.has_granular_temperature()) {

    const int ic_count = ics.size();

    std::vector<Real> meanVel(4*ic_count, 0);

    // First block: compute mean to force back to zero.
    for (int icv(0); icv < ics.size(); icv++) {

      // Avoid destroying fluctuations with nested regions
      if (!m_initial_conditions.has_granular_temperature(icv))
        continue;

      auto& ic = ics[icv];
      auto& ic_region = ic.region;

      const RealBox ic_realbox(ic_region->lo(), ic_region->hi());

      ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
      ReduceData<int, Real, Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
        Box tilebox = pti.tilebox();
        AoS& aos = pti.GetArrayOfStructs();
        ParticleType* pstruct = aos().dataPtr();

        RealBox tilebox_region(tilebox, Geom(lev).CellSize(), Geom(lev).ProbLo());

        if (tilebox_region.intersects(ic_realbox)) {

          const int np = pti.numParticles();

          SoA& soa = pti.GetStructOfArrays();
          auto p_realarray = soa.realarray();

          reduce_op.eval(np, reduce_data, [pstruct,p_realarray,ic_realbox]
            AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
          {
            const ParticleType p = pstruct[p_id];
            int cnt = 0;
            Real velx(0.);
            Real vely(0.);
            Real velz(0.);
            if ( ic_realbox.contains(p.pos()) ) {
              cnt = 1;
              velx = p_realarray[SoArealData::velx][p_id];
              vely = p_realarray[SoArealData::vely][p_id];
              velz = p_realarray[SoArealData::velz][p_id];
            }
            return {cnt, velx, vely, velz};
          });
        }
      }
      ReduceTuple host_tuple = reduce_data.value();

      meanVel[0*ic_count + icv] += static_cast<Real>(amrex::get<0>(host_tuple));

      meanVel[1*ic_count + icv] += amrex::get<1>(host_tuple);
      meanVel[2*ic_count + icv] += amrex::get<2>(host_tuple);
      meanVel[3*ic_count + icv] += amrex::get<3>(host_tuple);
    }

    ParallelDescriptor::ReduceRealSum(meanVel.data(), meanVel.size());

    std::vector<Real> granTemp(2*ic_count, 0.);

    // Subtract the mean velocities from each particle random velocity
    // so the new means are zero. Also, compute the mean granular
    // temperature.
    for (int icv(0); icv < ics.size(); icv++) {

      auto& ic = ics[icv];
      auto& ic_region = ic.region;

      const RealBox ic_realbox(ic_region->lo(), ic_region->hi());

      Real ic_np = meanVel[icv];

      if (ic_np > 0.) {

        const Real meanU = meanVel[1*ic_count + icv] / ic_np;
        const Real meanV = meanVel[2*ic_count + icv] / ic_np;
        const Real meanW = meanVel[3*ic_count + icv] / ic_np;


        ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<int, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

          Box tilebox = pti.tilebox();
          AoS& aos = pti.GetArrayOfStructs();
          ParticleType* pstruct = aos().dataPtr();

          RealBox tilebox_region(tilebox, Geom(lev).CellSize(), Geom(lev).ProbLo());

          if (tilebox_region.intersects(ic_realbox)) {

            const int np = pti.numParticles();

            SoA& soa = pti.GetStructOfArrays();
            auto p_realarray = soa.realarray();

            reduce_op.eval(np, reduce_data, [pstruct,p_realarray,ic_realbox,
              meanU, meanV, meanW]
              AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
            {
              const ParticleType p = pstruct[p_id];

              int cnt = 0;
              Real p_granTemp(0.);

              if ( ic_realbox.contains(p.pos()) ) {

                cnt = 1;

                Real velx = p_realarray[SoArealData::velx][p_id] - meanU;
                Real vely = p_realarray[SoArealData::vely][p_id] - meanV;
                Real velz = p_realarray[SoArealData::velz][p_id] - meanW;

                p_granTemp = velx*velx + vely*vely + velz*velz;

                p_realarray[SoArealData::velx][p_id] = velx;
                p_realarray[SoArealData::vely][p_id] = vely;
                p_realarray[SoArealData::velz][p_id] = velz;
              }
              return {cnt, p_granTemp};
            });
          } // tilebox intersects icv
        } // MFIter loop

          ReduceTuple host_tuple = reduce_data.value();

          granTemp[0*ic_count + icv] = static_cast<Real>(amrex::get<0>(host_tuple));
          granTemp[1*ic_count + icv] = amrex::get<1>(host_tuple);

      } // ic region has particles
    } // loop over ic regions

    ParallelDescriptor::ReduceRealSum(granTemp.data(), granTemp.size());

    for (int icv(0); icv < ics.size(); icv++) {

      auto& ic = ics[icv];
      auto& ic_region = ic.region;

      const RealBox ic_realbox(ic_region->lo(), ic_region->hi());

      Real ic_np = granTemp[icv];

      if (ic_np > 0.) {

        GpuArray<Real,3> bulkVel;

        auto& ic_solids = ic.solids;

        for (int lcs(0); lcs < ic_solids.size(); lcs++) {

          if (ic_solids[lcs].volfrac > tolerance) {

            const int phase = ic_solids[lcs].phase;

            const auto ic_solid = ic.get_solid(phase);

            if (ic_solid == nullptr) {
              amrex::Abort("Error");
            } else {
              bulkVel[0]= ic_solid->velocity[0];
              bulkVel[1]= ic_solid->velocity[1];
              bulkVel[2]= ic_solid->velocity[2];
            }

          }
        }

        const Real ic_granTemp(granTemp[ic_count + icv] / (3.0*ic_np));
        const Real scale(std::sqrt(m_initial_conditions.get_granular_temperature(icv)/ic_granTemp));

        // Adjust velocities so that the mean granular temperature is equal
        // to the desired initial granular temperature from the inputs.

        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

          Box tilebox = pti.tilebox();
          AoS& aos = pti.GetArrayOfStructs();
          ParticleType* pstruct = aos().dataPtr();

          RealBox tilebox_region(tilebox, Geom(lev).CellSize(), Geom(lev).ProbLo());

          if (tilebox_region.intersects(ic_realbox)) {

            const int np = pti.numParticles();

            SoA& soa = pti.GetStructOfArrays();
            auto p_realarray = soa.realarray();

            amrex::ParallelFor(np, [pstruct,p_realarray,ic_realbox, bulkVel, scale]
              AMREX_GPU_DEVICE (int p_id) noexcept
            {
              const ParticleType p = pstruct[p_id];

              if ( ic_realbox.contains(p.pos()) ) {

                Real velx = p_realarray[SoArealData::velx][p_id];
                Real vely = p_realarray[SoArealData::vely][p_id];
                Real velz = p_realarray[SoArealData::velz][p_id];

                p_realarray[SoArealData::velx][p_id] = bulkVel[0] + scale*velx;
                p_realarray[SoArealData::vely][p_id] = bulkVel[1] + scale*vely;
                p_realarray[SoArealData::velz][p_id] = bulkVel[2] + scale*velz;

             }
            });

          } // tilebox intersects ic region
        } // MFIter loop
      } // ic region has particles
    } // // loop over icv
  } // has granular temperature


}


void MFIXParticleContainer::
InitParticlesRuntimeVariables (const int adv_enthalpy)
{
  int lev = 0;
  const auto dx  = Geom(lev).CellSizeArray();
  const auto dx_inv = Geom(lev).InvCellSizeArray();

  const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();

  const int solve_species = m_solids.solve_species();

  const int nspecies_s = m_solids.nspecies();
  const int solid_is_a_mixture = m_solids.isMixture();

  const auto& solids_parms = m_solids.parameters<run_on>();

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

    auto& particles = pti.GetArrayOfStructs();
    int np = pti.numParticles();

    auto& soa = pti.GetStructOfArrays();
    auto p_realarray = soa.realarray();
    auto p_intarray = soa.intarray();

    auto particles_ptr = particles().dataPtr();

    PairIndex index(pti.index(), pti.LocalTileIndex());
    auto& plev = GetParticles(lev);
    auto& ptile = plev[index];
    auto ptile_data = ptile.getParticleTileData();

    // Initialize runtime reals that are not associated with an IC
    int const contains_eps = m_runtimeRealData.contains_eps();

    // ep_s
    if (contains_eps) {

      int const idx_eps = m_runtimeRealData.ep_s;

      ParallelFor(np, [ptile_data, idx_eps]
      AMREX_GPU_DEVICE (int ip) noexcept
      {
        ptile_data.m_runtime_rdata[idx_eps][ip] = 0.;
      });

    } // solve PIC

    int const contains_vm = m_runtimeRealData.contains_vm();
    int const contains_acc = m_runtimeRealData.contains_acc();

    if (contains_vm || contains_acc) {

      int const idx_vm_coeff = m_runtimeRealData.vm_coeff;
      int const idx_acceleration = m_runtimeRealData.acceleration;

      ParallelFor(np, [ptile_data, contains_vm, idx_vm_coeff,
        contains_acc, idx_acceleration]
      AMREX_GPU_DEVICE (int ip) noexcept
      {
        if ( contains_vm ) {
          ptile_data.m_runtime_rdata[idx_vm_coeff][ip] = 0.;
        }
        if ( contains_acc ) {
          ptile_data.m_runtime_rdata[idx_acceleration  ][ip] = 0.;
          ptile_data.m_runtime_rdata[idx_acceleration+1][ip] = 0.;
          ptile_data.m_runtime_rdata[idx_acceleration+2][ip] = 0.;
        }
      });

    } // contains vm or acc

    if (m_runtimeRealData.contains_energy()) {

      int const scomp = m_runtimeRealData.energy_source;
      int const ncomp = m_runtimeRealData.ncomp_energy;

      ParallelFor(np, [ptile_data, scomp, ncomp]
      AMREX_GPU_DEVICE (int ip) noexcept
      {
        for (int n(0); n<ncomp; ++n) {
          ptile_data.m_runtime_rdata[scomp + n][ip] = 0.;
        }
      });

    } // contains energy

    // Initialize tangential history data on all particles
    const int count_th_real = 3*(m_runtimeRealData.max_contacts);
    const int pft_neighbor_idx = m_runtimeRealData.pft_neighbor_idx;
    const int intcount = m_runtimeIntData.count;
    const int cpu_id_idx = m_runtimeIntData.cpu_id_idx;
    amrex::ParallelFor(np,
      [ptile_data,count_th_real,intcount,pft_neighbor_idx,cpu_id_idx]
      AMREX_GPU_DEVICE (int ip) noexcept
    {
        for (int n_th(0); n_th < count_th_real; ++n_th) {
           ptile_data.m_runtime_rdata[pft_neighbor_idx+n_th][ip] = 0.0;
        }

        for (int n_th(0); n_th < intcount; ++n_th) {
           ptile_data.m_runtime_idata[cpu_id_idx+n_th][ip] = -1;
        }
    });

    // Set the initial conditions.
    for(int icv(0); icv < m_initial_conditions.ic().size(); ++icv) {

      const IC_t& loc_ic = m_initial_conditions.ic(icv);

      // This is a round about way to address what is likely overly complex
      // logic in the particle generator. We take the region specified by
      // a user and convert it index space (i,j,k). That in turn is turned
      // back into a physical region that may be a little larger than what
      // was actually defined to account for spatial discretization.
      const IntVect bx_lo(static_cast<int>(amrex::Math::floor((loc_ic.region->lo(0)-plo[0])*dx_inv[0] + 0.5)),
                          static_cast<int>(amrex::Math::floor((loc_ic.region->lo(1)-plo[1])*dx_inv[1] + 0.5)),
                          static_cast<int>(amrex::Math::floor((loc_ic.region->lo(2)-plo[2])*dx_inv[0] + 0.5)));

      const IntVect bx_hi(static_cast<int>(amrex::Math::floor((loc_ic.region->hi(0)-plo[0])*dx_inv[0] + 0.5)),
                          static_cast<int>(amrex::Math::floor((loc_ic.region->hi(1)-plo[1])*dx_inv[1] + 0.5)),
                          static_cast<int>(amrex::Math::floor((loc_ic.region->hi(2)-plo[2])*dx_inv[0] + 0.5)));

      const Box ic_box(bx_lo, bx_hi);

      if (pti.tilebox().intersects(ic_box)) {

        // Start/end of IC domain bounds
        const RealVect ic_lo = {plo[0]+bx_lo[0]*dx[0],
                                plo[1]+bx_lo[1]*dx[1],
                                plo[2]+bx_lo[2]*dx[2]};

        const RealVect ic_hi = {plo[0]+bx_hi[0]*dx[0],
                                plo[1]+bx_hi[1]*dx[1],
                                plo[2]+bx_hi[2]*dx[2]};

        const RealBox ic_realbox(ic_lo.dataPtr(), ic_hi.dataPtr());

        // Loop through IC solids looking for match.
        for (int ics(0); ics < loc_ic.solids.size(); ics++) {

          auto& ic_solid = loc_ic.solids[ics];

          const int ic_phase = ic_solid.phase;

          // Create a temporary copy of IC particle temperatures mapped
          // to the particle type.
          Real h_temperature_loc(0.);
          Gpu::HostVector<Real> h_mass_fractions(nspecies_s);

          if (adv_enthalpy) {
            h_temperature_loc = ic_solid.temperature;
          }

          if (solve_species) {
            for (int n_s(0); n_s < nspecies_s; n_s++)
              h_mass_fractions[n_s] = ic_solid.species[n_s].mass_fraction;
          }

          Gpu::AsyncArray<Real> d_mass_fractions(h_mass_fractions.data(), h_mass_fractions.size());
          Real* p_mass_fractions = solve_species ? d_mass_fractions.data() : nullptr;

          const int idx_X_sn = m_runtimeRealData.X_sn;

          amrex::ParallelFor(np,
            [particles_ptr,p_realarray,p_intarray,ptile_data,h_temperature_loc,
             p_mass_fractions,ic_realbox,nspecies_s,solid_is_a_mixture,adv_enthalpy,
             solids_parms,solve_species,idx_X_sn,ic_phase]
            AMREX_GPU_DEVICE (int ip) noexcept
          {
            const auto& p = particles_ptr[ip];

            const int p_phase = p_intarray[SoAintData::phase][ip];

            if(ic_realbox.contains(p.pos()) && (p_phase == ic_phase)) {

              if(adv_enthalpy) {
                p_realarray[SoArealData::temperature][ip] = h_temperature_loc;
              }

              if(solve_species) {
                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  ptile_data.m_runtime_rdata[idx_X_sn+n_s][ip] = p_mass_fractions[n_s];
                }
              }
            }
          });

        } // for ic_solids.size()


      } // Intersecting Boxes
    } // IC regions
  } // MFIXParIter
}
