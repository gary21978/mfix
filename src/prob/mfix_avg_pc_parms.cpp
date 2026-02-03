#include <mfix_avg_pc_parms.H>
#include <mfix_solids.H>

#include <mfix_reporter.H>

using namespace amrex;


void MFIXAvgParticleParms::
Initialize ( int const a_pc_var_count,
             MFIXFluidPhase& a_fluid,
             MFIXEmbeddedBoundaries& a_eb_parms)
{
  m_pc_var_count = a_pc_var_count;

  m_h_depComp.resize(a_pc_var_count, 0);
  m_comp_map.resize(a_pc_var_count, -1);

  if ( a_fluid.susp_visc_model( SuspViscModel::Sato ) ||
       a_fluid.susp_visc_model( SuspViscModel::Subramaniam ) ) {

    setDepComp(SoArealData::radius);
    setDepComp(SoArealData::velx);
    setDepComp(SoArealData::vely);
    setDepComp(SoArealData::velz);
  }

  if ( a_eb_parms.needs_avg_particle_temperature() ) {
    setDepComp(SoArealData::temperature);
  }
}
