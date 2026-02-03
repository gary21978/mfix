#include <AMReX_ParmParse.H>

#include <pm_particles.H>
#include <pm_reporter.H>

using namespace amrex;

pm_particles::
pm_particles ( int const a_max_level,
               Vector<Geometry>            const & a_geom,
               Vector<DistributionMapping> const & a_dmap,
               Vector<BoxArray>            const & a_ba,
               pm_plotfile* a_plotfile)
  : amrex::NeighborParticleContainer<0, 0, 0, 0>(a_geom[0], a_dmap[0], a_ba[0], 1)
  , m_max_level(a_max_level)
  , m_verbose(0)
  , m_real_comps(a_plotfile->particle_real_comps())
  , m_int_comps(a_plotfile->particle_int_comps())
  , m_int_comp_names(a_plotfile->particle_int_variables())
  , m_real_comp_names(a_plotfile->particle_real_variables())
  , m_mean_diameter(-1.0)
{

  this->SetVerbose(0);

  for( int n(0); n<m_int_comps; ++n) {
    AddIntComp(true);
    m_variable_map_i[m_int_comp_names[n]] = n;
  }

  for( int n(0); n<m_real_comps; ++n) {
    AddRealComp(true);
    m_variable_map_r[m_real_comp_names[n]] = n;
  }

  // Positions
  setRealCommComp(0, true);
  setRealCommComp(1, true);
  setRealCommComp(2, true);

  Restart(a_plotfile->get_name(), "particles");

  if ( contains("radius") ) {
    int const radius_index = var_index("radius") + AMREX_SPACEDIM;
    setRealCommComp(radius_index, true);
  }

  { ParmParse pp_particle("particle");

    if ( pp_particle.query("mean_diameter", m_mean_diameter) ) {
      Print() << "User provided mean diameter: " << m_mean_diameter << "\n";

    } else if ( contains("radius") ) {
      AMREX_ALWAYS_ASSERT( contains("radius") );

      Real inv_pcount = 1.0/static_cast<Real>(TotalNumberOfParticles());
      m_mean_diameter = (2.0*sum("radius")) * inv_pcount;

      Print() << "Computed mean diameter: " << m_mean_diameter << "\n";
    }

  }

  if ( !contains("volume") && !contains("radius") ) {
    post_mfix::Error(__FILE__, __LINE__) << " Error:"
      << " Plotfile does not contain particle radius or volume.\n"
      << " Unable to proceed.";
  }

}



Real pm_particles::
max ( const std::string a_name,
      bool const a_local )
{
  int const comp = var_index(a_name);

  Real mx(-1.);

  for (int lev(0); lev < get_nlev(); lev++) {

    auto& plev = GetParticles(lev);

    // Reduce sum operation
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (PMParIter pti(*this, lev); pti.isValid(); ++pti) {

      std::pair<int, int> index(pti.index(), pti.LocalTileIndex());

      // variables added at runtime
      auto ptile_data = plev[index].getParticleTileData();

      reduce_op.eval(pti.numParticles(), reduce_data, [ptile_data, comp]
      AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
      { return {ptile_data.m_runtime_rdata[comp][p_id]}; });
    } // MFIXParIter

    ReduceTuple host_tuple = reduce_data.value();
    mx = amrex::max(mx, amrex::get<0>(host_tuple));

  }

  if (!a_local) { ParallelAllReduce::Max(mx, ParallelContext::CommunicatorSub()); }

  return mx;
}


Real pm_particles::
sum ( const std::string a_name,
      bool const a_local )
{
  int const comp = var_index(a_name);

  Real sm(0.);

  for (int lev(0); lev < get_nlev(); lev++) {

    auto& plev = GetParticles(lev);

    // Reduce sum operation
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (PMParIter pti(*this, lev); pti.isValid(); ++pti) {

      std::pair<int, int> index(pti.index(), pti.LocalTileIndex());

      // variables added at runtime
      auto ptile_data = plev[index].getParticleTileData();

      reduce_op.eval(pti.numParticles(), reduce_data, [ptile_data, comp]
      AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
      { return {ptile_data.m_runtime_rdata[comp][p_id]}; });
    } // MFIXParIter

    ReduceTuple host_tuple = reduce_data.value();
    sm += amrex::get<0>(host_tuple);

  }

  if (!a_local) { ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub()); }

  return sm;
}


Real pm_particles::
fluct ( std::string const a_var,
        Real const a_avg,
        bool const a_local )
{
  int const comp = var_index(a_var);

  Real sm(0.);

  for (int lev(0); lev < get_nlev(); lev++) {

    auto& plev = GetParticles(lev);

    // Reduce sum operation
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (PMParIter pti(*this, lev); pti.isValid(); ++pti) {

      std::pair<int, int> index(pti.index(), pti.LocalTileIndex());

      // variables added at runtime
      auto ptile_data = plev[index].getParticleTileData();

      reduce_op.eval(pti.numParticles(), reduce_data, [ptile_data, comp, a_avg]
      AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
      { Real t = ptile_data.m_runtime_rdata[comp][p_id] - a_avg;
        return t*t;
      });
    } // MFIXParIter

    ReduceTuple host_tuple = reduce_data.value();
    sm += amrex::get<0>(host_tuple);

  }

  if (!a_local) { ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub()); }

  return sm;
}
