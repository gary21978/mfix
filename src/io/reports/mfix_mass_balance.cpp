#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>

#include <mfix_rw.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_monitors.H>
#include <mfix_reporter.H>
#include <mfix_mass_balance.H>

using namespace amrex;


MFIXMassBalance::
MFIXMassBalance (amrex::Vector<amrex::Geometry> const& a_geom,
                 amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& a_ebfactory,
                 MFIXLevelData& a_leveldata,
                 MFIXParticleContainer* a_pc,
                 MFIXFluidPhase const& a_fluid,
                 MFIXSolidsPhase const& a_solids,
                 MFIXDEM const& a_dem,
                 MFIXPIC const& a_pic,
                 BCList const& a_bclist)
  : m_nlev(a_geom.size())
  , m_geom(a_geom)
  , m_ebfactory(a_ebfactory)
  , m_leveldata(a_leveldata)
  , m_pc(a_pc)
  , m_fluid(a_fluid)
  , m_solids(a_solids)
  , m_dem(a_dem)
  , m_pic(a_pic)
  , m_bclist(a_bclist)
{}


void
MFIXMassBalance::
Initialize ()
{
  m_fluid_mass_accum.resize(m_fluid.nspecies()*2, 0.);
  m_fluid_mass_inflow.resize(m_fluid.nspecies(), 0.);
  m_fluid_mass_outflow.resize(m_fluid.nspecies(), 0.);
  m_fluid_mass_prod.resize(m_fluid.nspecies(), 0.);

  m_solids_mass_accum.resize(m_solids.nspecies()*2, 0.);
  m_solids_mass_inflow.resize(m_solids.nspecies(), 0.);
  m_solids_mass_outflow.resize(m_solids.nspecies(), 0.);
  m_solids_mass_prod.resize(m_solids.nspecies(), 0.);
}


void
MFIXMassBalance::
ComputeMassAccum (const int report_mass_balance,
                  const int a_offset)
{
  BL_PROFILE("mfix::ComputeMassAccum()");

  if (!report_mass_balance) {
    return;
  }

  if (m_fluid.solve()) {

    Vector<Real> accum(m_fluid.nspecies(), 0.);

    const int nspecies_g = m_fluid.nspecies();
    for (int lev = 0; lev < m_nlev; lev++) {

      const GpuArray<Real,3> dx = m_geom[lev].CellSizeArray();
      const Real vol = dx[0]*dx[1]*dx[2];

      for (MFIter mfi(*(leveldata().X(lev)), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        EBCellFlagFab const& flagfab = m_ebfactory[lev]->getMultiEBCellFlagFab()[mfi];

        if (flagfab.getType(bx) != FabType::covered) {

          Array4<Real const> const& epf = leveldata().epf_const(lev,mfi);
          Array4<Real const> const& rho = leveldata().rho_const(lev,mfi);
          Array4<Real const> const& X_n = leveldata().X_const(lev,mfi);

          if (flagfab.getType(bx) == FabType::singlevalued ) {

            Array4<Real const> const& volfrac = (m_ebfactory[lev]->getVolFrac()).const_array(mfi);

            for (int n=0; n < nspecies_g; ++n){

              ReduceOps<ReduceOpSum> reduce_op;
              ReduceData<Real> reduce_data(reduce_op);
              using ReduceTuple = typename decltype(reduce_data)::Type;

              reduce_op.eval(bx, reduce_data, [n, epf, rho, X_n, volfrac]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                return { volfrac(i,j,k)*epf(i,j,k)*rho(i,j,k)*X_n(i,j,k,n) };
              });

              ReduceTuple host_tuple = reduce_data.value(reduce_op);
              accum[n] += amrex::get<0>(host_tuple)*vol;

            } /* End loop over species */

          } else {

            for (int n=0; n < nspecies_g; ++n){

              ReduceOps<ReduceOpSum> reduce_op;
              ReduceData<Real> reduce_data(reduce_op);
              using ReduceTuple = typename decltype(reduce_data)::Type;

              reduce_op.eval(bx, reduce_data, [n, epf, rho, X_n]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                return { epf(i,j,k)*rho(i,j,k)*X_n(i,j,k,n) };
              });

              ReduceTuple host_tuple = reduce_data.value(reduce_op);
              accum[n] += amrex::get<0>(host_tuple)*vol;

            } /* End loop over species */
          } // FabType
        } // Not covered Fab
      } // loop over MFIter

    } // m_nlev

    // Global sum and copy to global variable with a_offset
    ParallelDescriptor::ReduceRealSum(accum.dataPtr(), nspecies_g);
    for (int n=0; n < nspecies_g; ++n) {
      m_fluid_mass_accum[n + a_offset*nspecies_g] = accum[n];
    }

  } // solve fluid

  if (m_dem.solve() || m_pic.solve()) {
    using MFIXParIter = MFIXParticleContainer::MFIXParIter;
    using PairIndex = MFIXParticleContainer::PairIndex;

    const int nspecies_s = m_solids.nspecies();
    std::vector<Real> accum(nspecies_s, 0.);

    for (int lev(0); lev < m_nlev; lev++) {

      const int idx_X_sn = m_pc->m_runtimeRealData.X_sn;

      // Ideally, we would do this in one loop, but since the number of species
      // is unknown, we loop over each one.
      for (int n(0); n < nspecies_s; ++n) {

        // Reduce sum operation for mass and production of n-th species
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIXParIter pti(*m_pc, lev); pti.isValid(); ++pti) {

          auto& plev  = m_pc->GetParticles(lev);
          PairIndex index(pti.index(), pti.LocalTileIndex());

          // SoA real variables
          auto& soa = pti.GetStructOfArrays();

          // particle stat weight and mass
          auto p_statwt = soa.GetRealData(SoArealData::statwt).data();
          auto p_radius = soa.GetRealData(SoArealData::radius).data();
          auto p_density = soa.GetRealData(SoArealData::density).data();

          // variables added at runtime
          auto ptile_data = plev[index].getParticleTileData();

          const int np = pti.numParticles();

          reduce_op.eval(np, reduce_data, [p_statwt, p_radius, p_density,
              ptile_data, idx_X_sn, n]
          AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
          {
            const Real p_mass = p_density[p_id]*SoArealData::volume(p_radius[p_id]);
            return {p_statwt[p_id]*p_mass*ptile_data.m_runtime_rdata[idx_X_sn+n][p_id]};
          });
        } // MFIXParIter

        ReduceTuple host_tuple = reduce_data.value();
        accum[n] = amrex::get<0>(host_tuple);

      } //loop over species

    }//lev

    ParallelDescriptor::ReduceRealSum(accum.data(), nspecies_s);
    for (int n=0; n < nspecies_s; ++n) {
      m_solids_mass_accum[n + a_offset*nspecies_s] += accum[n];
    }
  }

}


void
MFIXMassBalance::ComputeMassProduction (Gpu::DeviceVector<Real> const& d_fluid_prod,
                                        Gpu::DeviceVector<Real> const& d_solids_prod)
{
  const int nspecies_g = m_fluid.nspecies();
  const int nspecies_s = m_solids.nspecies();

  // Global sum and copy to global variable
  Gpu::HostVector<Real> h_fluid_prod(nspecies_g, 0.);
  Gpu::HostVector<Real> h_solids_prod(nspecies_s, 0.);

#ifdef AMREX_USE_GPU
  Gpu::copy(Gpu::deviceToHost, d_fluid_prod.begin(), d_fluid_prod.end(), h_fluid_prod.begin());
  Gpu::copy(Gpu::deviceToHost, d_solids_prod.begin(), d_solids_prod.end(), h_solids_prod.begin());
#else
  std::copy(d_fluid_prod.begin(), d_fluid_prod.end(), h_fluid_prod.begin());
  std::copy(d_solids_prod.begin(), d_solids_prod.end(), h_solids_prod.begin());
#endif

  ParallelDescriptor::ReduceRealSum(h_fluid_prod.data(), nspecies_g);
  ParallelDescriptor::ReduceRealSum(h_solids_prod.data(), nspecies_s);

  // Global sum and copy to global variable
  Gpu::HostVector<Real>& fluid_mass_prod = this->fluid_mass_prod();
  Gpu::HostVector<Real>& solids_mass_prod = this->solids_mass_prod();

  std::transform(fluid_mass_prod.begin(), fluid_mass_prod.end(), h_fluid_prod.begin(),
      fluid_mass_prod.begin(), std::plus<Real>());

  std::transform(solids_mass_prod.begin(), solids_mass_prod.end(), h_solids_prod.begin(),
      solids_mass_prod.begin(), std::plus<Real>());
}


void
MFIXMassBalance::ComputeMassFlux (Vector< MultiFab const*> const& flux_x,
                                  Vector< MultiFab const*> const& flux_y,
                                  Vector< MultiFab const*> const& flux_z,
                                  const int scomp,
                                  const int ncomp,
                                  const bool fluxes_are_area_weighted,
                                  const int eb_has_flow,
                                  Vector< MultiFab const*> const& eb_vel_in,
                                  Vector< MultiFab const*> const& eb_species_in,
                                  const Real dt)
{
  amrex::ignore_unused(ncomp, fluxes_are_area_weighted);

  const int nspecies_g = m_fluid.nspecies();

  std::vector<Real> mass_flow(2*nspecies_g, 0.);

  for (int lev = 0; lev < m_nlev; lev++) {

    Array<const MultiFab*,3> flux = {flux_x[lev], flux_y[lev], flux_z[lev]};

    const GpuArray<Real,3> dx = m_geom[lev].CellSizeArray();

    // This should be caught elsewhere but just in case...
    AMREX_ASSERT(dx[0] == dx[1] && dx[1] == dx[2]);
    Real const da( dx[0]*dx[0] );

    amrex::EBFArrayBoxFactory const& fact =
      static_cast<amrex::EBFArrayBoxFactory const&>(*m_ebfactory[lev]);

    Box domain(m_geom[lev].Domain());

    Array<Array4<int>,3> bct_lo = {m_bclist.bc_ilo[lev]->array(),
                                   m_bclist.bc_jlo[lev]->array(),
                                   m_bclist.bc_klo[lev]->array()};

    Array<Array4<int>,3> bct_hi = {m_bclist.bc_ihi[lev]->array(),
                                   m_bclist.bc_jhi[lev]->array(),
                                   m_bclist.bc_khi[lev]->array()};

    Array<Array<Array4<int>,3>,2> bct_extents = {bct_lo, bct_hi};

    auto areafrac = fact.getAreaFrac();

    for (int n=0; n < nspecies_g; ++n) {

      ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
      ReduceData<Real,Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      // Loop over MFIter
      for (MFIter mfi(*(leveldata().epf(lev)), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];

        Array<IntVect,2> domain_extents = {domain.smallEnd(), domain.bigEnd()};

        // Loop over directions
        for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

          // Get flux over direction dir
          Array4<Real const> const& flux_arr = flux[dir]->const_array(mfi);

          const Box& sbox = (*flux[dir])[mfi].box();

          // Save local box extents in a specific container
          Array<IntVect,2> sbox_extents = {sbox.smallEnd(), sbox.bigEnd()};

          // Loop over face normal
          for (int normal(-1); normal <= 1; normal += 2) {

            // index to map normal from {-1,1} to {0,1}
            const int idx = (normal + 1) / 2;

            // On domain face
            if (sbox_extents[idx][dir] == (domain_extents[idx][dir]+idx)) {

              // Create a copy of local box extents
              Array<IntVect,2> operative_box_extents(sbox_extents);

              // Modify operative box extents to make a 2D Box
              operative_box_extents[(idx+1)%2][dir] = operative_box_extents[idx][dir];

              // Create the 2D Box
              const Box operative_box(operative_box_extents[0],
                                      operative_box_extents[1]);

              // Get the boundary condition type on current face
              Array4<int>& bc_type = bct_extents[idx][dir];

              // Get domain index on this face
              const int domain_idx = domain_extents[idx][dir];

              if (flagfab.getType(bx) == FabType::regular ) {

                reduce_op.eval(operative_box, reduce_data, [dir,bc_type,
                domain_idx,normal,flux_arr,n,scomp]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                  Real m_in  = 0.;
                  Real m_out = 0.;

                  IntVect bct_cell(i,j,k);
                  bct_cell[dir] = domain_idx + normal;
                  const int bct = bc_type(bct_cell,0);

                  if(bct == BCList::pout) {
                    //m_out =  1.;
                    m_out =  normal*flux_arr(i,j,k,n+scomp);
                  }
                  else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                    //m_in  = 1.;
                    m_in  = -normal*flux_arr(i,j,k,n+scomp);
                  }

                  return {m_in, m_out};
                });


              } else if (flagfab.getType(bx) == FabType::singlevalued ) {

                Array4<Real const> const& areafrac_arr = areafrac[dir]->const_array(mfi);

                reduce_op.eval(operative_box, reduce_data, [dir,bc_type,
                domain_idx,normal,flux_arr,areafrac_arr,n,scomp]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                  Real m_in  = 0.;
                  Real m_out = 0.;

                  IntVect bct_cell(i,j,k);
                  bct_cell[dir] = domain_idx + normal;
                  const int bct = bc_type(bct_cell,0);

                  if(bct == BCList::pout) {
                    m_out =  normal*areafrac_arr(i,j,k)*flux_arr(i,j,k,n+scomp);
                    //m_out =  areafrac_arr(i,j,k);
                  }
                  else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                    m_in  = -normal*areafrac_arr(i,j,k)*flux_arr(i,j,k,n+scomp);
                    //m_in  = areafrac_arr(i,j,k);
                  }

                  return {m_in, m_out};
                });

              } // single valued fab
            } // on domain face
          } // loop over normal direction
        } // loop over direction

        // Add in mass flow from EB
        if( eb_has_flow && flagfab.getType(bx) == FabType::singlevalued) {

          // Area of eb face
          Array4<Real const> const& barea  = fact.getBndryArea().const_array(mfi);

          Array4<Real const> const& eb_vel     = eb_vel_in[lev]->const_array(mfi);
          Array4<Real const> const& eb_species = eb_species_in[lev]->const_array(mfi);

          reduce_op.eval(bx, reduce_data, [n, eb_vel, eb_species, barea]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
          {

            Real vel_mag = std::sqrt(eb_vel(i,j,k,0)*eb_vel(i,j,k,0) +
                                     eb_vel(i,j,k,1)*eb_vel(i,j,k,1) +
                                     eb_vel(i,j,k,2)*eb_vel(i,j,k,2));

            return {vel_mag*eb_species(i,j,k,n)*barea(i,j,k), 0.};

          });

        } // eb_has_flow

      } // MFIter loop

      ReduceTuple host_tuple = reduce_data.value(reduce_op);

      // mass in
      mass_flow[n] += amrex::get<0>(host_tuple)*da;
      // mass out
      mass_flow[n+nspecies_g] += amrex::get<1>(host_tuple)*da;

    } // Loop over species

  } // loop over levels

  ParallelDescriptor::ReduceRealSum(mass_flow.data(), 2*nspecies_g);

  // Copy into global variables.
  for (int n=0; n < nspecies_g; ++n) {

    m_fluid_mass_inflow[n]  += dt*mass_flow[n];
    m_fluid_mass_outflow[n] += dt*mass_flow[n+nspecies_g];
  }
}


void
MFIXMassBalance::InitMassBalance (const int report_mass_balance,
                                  const int nspecies)
{
  if (report_mass_balance) {
    for(int n(0); n < nspecies; n++) {
        m_fluid_mass_accum[n] = m_fluid_mass_accum[n+nspecies];
        m_fluid_mass_inflow[n] = 0.;
        m_fluid_mass_outflow[n] = 0.;
        m_fluid_mass_prod[n] = 0.;
    }
  }
}


void
MFIXMassBalance::ResetMassBalance (int const a_n)
{
  m_solids_mass_accum[a_n] = m_solids_mass_accum[a_n + m_solids.nspecies()];

  m_solids_mass_accum[a_n + m_solids.nspecies()] = 0.;
  m_solids_mass_inflow[a_n]  = 0.;
  m_solids_mass_outflow[a_n] = 0.;
  m_solids_mass_prod[a_n]    = 0.;
}


void
MFIXMassBalance::ComputeMassOutflow (int const a_lev)
{
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;
  using PairIndex = MFIXParticleContainer::PairIndex;

  const int nspecies_s = m_solids.nspecies();
  std::vector<Real> outflow(nspecies_s, 0.);

  // particle tiles and geometry of this level

  const auto p_lo = m_geom[a_lev].ProbLoArray();
  const auto p_hi = m_geom[a_lev].ProbHiArray();

  const int idx_X_sn = m_pc->m_runtimeRealData.X_sn;

  // Ideally, we would do this in one loop, but since the number of species
  // is unknown, we loop over each one.
  for (int n(0); n < nspecies_s; ++n) {

    // Reduce sum operation for mass and production of n-th species
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIXParIter pti(*m_pc, a_lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& plev  = m_pc->GetParticles(a_lev);

      auto& aos   = plev[index].GetArrayOfStructs();
      MFIXParticleContainer::ParticleType* pstruct = aos().dataPtr();

      // SoA real variables
      auto& soa = pti.GetStructOfArrays();

      // particle stat weight and mass
      auto p_statwt = soa.GetRealData(SoArealData::statwt).data();
      auto p_radius   = soa.GetRealData(SoArealData::radius).data();
      auto p_density  = soa.GetRealData(SoArealData::density).data();

      // variables added at runtime
      auto ptile_data = plev[index].getParticleTileData();

      const int np = pti.numParticles();

      reduce_op.eval(np, reduce_data, [p_lo, p_hi, pstruct, p_statwt, p_radius,
          p_density, ptile_data, idx_X_sn, n]
      AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
      {
        auto p = pstruct[p_id];
        const Real p_mass = p_density[p_id]*SoArealData::volume(p_radius[p_id]);

        Real p_mass_Xn(0);
        if ( p.pos(0) < p_lo[0] || p.pos(0) > p_hi[0] ||
             p.pos(1) < p_lo[1] || p.pos(1) > p_hi[1] ||
             p.pos(2) < p_lo[2] || p.pos(2) > p_hi[2] ) {
          p_mass_Xn = p_statwt[p_id]*p_mass*ptile_data.m_runtime_rdata[idx_X_sn+n][p_id];
        }
        return {p_mass_Xn};
      });
    } // MFIXParIter

  ReduceTuple host_tuple = reduce_data.value();
  outflow[n] = amrex::get<0>(host_tuple);

  } //loop over species

  // Global sum and copy to global variable
  ParallelDescriptor::ReduceRealSum(outflow.data(), nspecies_s);
  for (int n=0; n < nspecies_s; ++n) {
    m_solids_mass_outflow[n] += outflow[n];
  }

}
