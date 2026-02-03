#include <AMReX_ParticleMesh.H>
#include <AMReX_ParticleInterpolators.H>

#include <pm_particles.H>

using namespace amrex;

void pm_particles::
deposit ( int const a_lev,
          MultiFab* a_MF, int const a_dstcomp,
          std::string const a_var)
{
  int var_comp = -1;

  if ( contains(a_var) ) { var_comp = var_index(a_var); }
  else { AMREX_ASSERT( toLower(a_var).compare("volume") == 0); }

  int const volume_comp = contains("volume") ? var_index("volume") : -1;
  int const radius_comp = contains("radius") ? var_index("radius") :  0;

  AMREX_ASSERT( contains("volume") || contains("radius"));

  a_MF->setVal(0.0, a_dstcomp, 1, a_MF->nGrow());

  const auto plo = Geom(a_lev).ProbLoArray();
  const auto dxi = Geom(a_lev).InvCellSizeArray();

  using CPTD_Type = pm_particles::ParticleTileType::ConstParticleTileDataType;

  amrex::ParticleToMesh(*this, *a_MF, a_lev,
  [plo, dxi, volume_comp, radius_comp, var_comp, a_dstcomp]
  AMREX_GPU_DEVICE (const CPTD_Type& ptd, int i, Array4<Real> const& MF)
  {

    // particle radius -- 0th component if radius not available in which case
    // we expect to use the stored particle volume.
    Real const radius = ptd.m_runtime_rdata[radius_comp][i];

    const amrex::Particle<0, 0> p = ptd.m_aos[i];
    ParticleInterpolator::Linear interp(p, plo, dxi);

    // Get the component value
    Real value = (var_comp < 0) ? 1.0 : ptd.m_runtime_rdata[var_comp][i];

    // Multiply component value by particle volume.
    value *= (volume_comp < 0) ? (4.0/3.0)*M_PI*(radius*radius*radius):
                                 ptd.m_runtime_rdata[volume_comp][i];

    int const start_particle_comp(0);
    int const start_mesh_comp(a_dstcomp);
    int const num_comp(1);

    interp.ParticleToMesh(p, MF, start_particle_comp, start_mesh_comp, num_comp,
    [=] AMREX_GPU_DEVICE (const pm_particles::ParticleType& /*part*/, int /*comp*/)
    { return value; });

  }, false); // false to prevent zeroing other components

  Real const inv_cell_vol = dxi[0]*dxi[1]*dxi[2];
  a_MF->mult(inv_cell_vol, a_dstcomp, 1, 0);
}


// Deposit x*y
void pm_particles::
deposit_mult ( int const a_lev,
               MultiFab* a_MF, int const a_dstcomp,
               std::string const a_var1,
               std::string const a_var2)
{
  AMREX_ASSERT( contains(a_var1) );
  int const var1_comp = var_index(a_var2);

  AMREX_ASSERT( contains(a_var2) );
  int const var2_comp = var_index(a_var2);

  int const volume_comp = contains("volume") ? var_index("volume") : -1;
  int const radius_comp = contains("radius") ? var_index("radius") :  0;

  AMREX_ASSERT( contains("volume") || contains("radius") );

  a_MF->setVal(0.0, a_dstcomp, 1, a_MF->nGrow());

  const auto plo = Geom(a_lev).ProbLoArray();
  const auto dxi = Geom(a_lev).InvCellSizeArray();

  using CPTD_Type = pm_particles::ParticleTileType::ConstParticleTileDataType;

  amrex::ParticleToMesh(*this, *a_MF, a_lev,
  [plo, dxi, volume_comp, radius_comp, var1_comp, var2_comp, a_dstcomp]
  AMREX_GPU_DEVICE (const CPTD_Type& ptd, int i, Array4<Real> const& MF)
  {
    const amrex::Particle<0, 0> p = ptd.m_aos[i];
    ParticleInterpolator::Linear interp(p, plo, dxi);

    Real const radius = ptd.m_runtime_rdata[radius_comp][i];

    // Multiply the components
    Real value = ptd.m_runtime_rdata[var1_comp][i] *
                 ptd.m_runtime_rdata[var2_comp][i];

    // Multiply component value by particle volume.
    value *= (volume_comp < 0) ? (4.0/3.0)*M_PI*(radius*radius*radius):
                                 ptd.m_runtime_rdata[volume_comp][i];

    int const start_particle_comp(0);
    int const start_mesh_comp(a_dstcomp);
    int const num_comp(1);

    interp.ParticleToMesh(p, MF, start_particle_comp, start_mesh_comp, num_comp,
    [=] AMREX_GPU_DEVICE (const pm_particles::ParticleType& /*part*/, int /*comp*/)
    { return value; });

  }, false); // false to prevent zeroing other components

  Real const inv_cell_vol = dxi[0]*dxi[1]*dxi[2];
  a_MF->mult(inv_cell_vol, a_dstcomp, 1, 0);

}


// Deposit a*(x^b)
void pm_particles::
deposit_pow ( int const a_lev,
              MultiFab* a_MF, int const a_dstcomp,
              std::string const a_var,
              Real const a_a, int const a_b)
{
  AMREX_ASSERT( contains(a_var) );
  int const var_comp = var_index(a_var);

  int const volume_comp = contains("volume") ? var_index("volume") : -1;
  int const radius_comp = contains("radius") ? var_index("radius") :  0;

  AMREX_ASSERT( contains("volume") || contains("radius") );

  a_MF->setVal(0.0, a_dstcomp, 1, a_MF->nGrow());

  const auto plo = Geom(a_lev).ProbLoArray();
  const auto dxi = Geom(a_lev).InvCellSizeArray();

  using CPTD_Type = pm_particles::ParticleTileType::ConstParticleTileDataType;

  amrex::ParticleToMesh(*this, *a_MF, a_lev,
  [plo, dxi, volume_comp, radius_comp, var_comp, a_a, a_b, a_dstcomp]
  AMREX_GPU_DEVICE (const CPTD_Type& ptd, int i, Array4<Real> const& MF)
  {
    const amrex::Particle<0, 0> p = ptd.m_aos[i];
    ParticleInterpolator::Linear interp(p, plo, dxi);

    Real const radius = ptd.m_runtime_rdata[radius_comp][i];

    // Evaluate the power: a*x^b
    Real value = a_a * std::pow(ptd.m_runtime_rdata[var_comp][i], a_b);

    // Multiply component value by particle volume.
    value *= (volume_comp < 0) ? (4.0/3.0)*M_PI*(radius*radius*radius):
                                 ptd.m_runtime_rdata[volume_comp][i];

    int const start_particle_comp(0);
    int const start_mesh_comp(a_dstcomp);
    int const num_comp(1);

    interp.ParticleToMesh(p, MF, start_particle_comp, start_mesh_comp, num_comp,
    [=] AMREX_GPU_DEVICE (const pm_particles::ParticleType& /*part*/, int /*comp*/)
    { return value; });

  }, false); // false to prevent zeroing other components

  Real const inv_cell_vol = dxi[0]*dxi[1]*dxi[2];
  a_MF->mult(inv_cell_vol, a_dstcomp, 1, 0);

}
