#include <mfix_stepdata.H>

using namespace amrex;

UniqueStepData::
UniqueStepData ( BoxArray const& a_ba, DistributionMapping const& a_dm,
                 FabFactory<FArrayBox> const& a_eb,
                 int const a_enthalpy,
                 int const a_nspecies,
                 int const a_ntracers,
                 int const a_ncomp_adv )
{
  m_conv_velocity.reset( new MultiFab(a_ba, a_dm, 3, 0, MFInfo(), a_eb) );
  m_conv_velocity->setVal(0.0);

  m_conv_scalars.reset( new MultiFab(a_ba, a_dm, a_ncomp_adv, 0, MFInfo(), a_eb) );
  m_conv_scalars->setVal(0.0);

  if ( a_enthalpy ) {
    m_lap_T.reset( new MultiFab(a_ba, a_dm, 1, 0, MFInfo(), a_eb) );
    m_lap_T->setVal(0.0);
  }

  if ( a_ntracers > 0 ) {
    m_lap_tracer.reset( new MultiFab(a_ba, a_dm, a_ntracers, 0, MFInfo(), a_eb) );
    m_lap_tracer->setVal(0.0);
  }

  if ( a_nspecies > 0 ) {
    m_conv_species.reset( new MultiFab(a_ba, a_dm, a_nspecies, 0, MFInfo(), a_eb) );
    m_conv_species->setVal(0.0);

    m_div_J.reset( new MultiFab(a_ba, a_dm, a_nspecies, 0, MFInfo(), a_eb) );
    m_div_J->setVal(0.0);
  }

  if ( a_enthalpy && (a_nspecies > 0) ) {
    m_div_hJ.reset( new MultiFab(a_ba, a_dm, 1, 0, MFInfo(), a_eb) );
    m_div_hJ->setVal(0.0);
  }

}

SharedStepData::
SharedStepData ( BoxArray const& a_ba, DistributionMapping const& a_dm,
                 FabFactory<FArrayBox> const& a_eb,
                 int const a_enthalpy,
                 int const a_nspecies,
                 int const a_ntracers,
                 int const /*a_ncomp_adv*/,
                 int const a_nghost_state,
                 int const a_nghost_force,
                 int const a_nghost_mac)
{
  const BoxArray& u_ba = amrex::convert(a_ba, IntVect::TheDimensionVector(0));
  const BoxArray& v_ba = amrex::convert(a_ba, IntVect::TheDimensionVector(1));
  const BoxArray& w_ba = amrex::convert(a_ba, IntVect::TheDimensionVector(2));

  m_vel_S_p.reset( new MultiFab(a_ba, a_dm, 1, 1, MFInfo(), a_eb) );
  m_vel_S_p->setVal(0.0);

  m_vel_S_c.reset( new MultiFab(a_ba, a_dm, 3, 1, MFInfo(), a_eb) );
  m_vel_S_c->setVal(0.0);

  m_vel_forces.reset( new MultiFab(a_ba, a_dm, 3, a_nghost_force, MFInfo(), a_eb) );
  m_vel_forces->setVal(0.0);

  if (a_ntracers > 0) {
    m_tracer_forces.reset(
        new MultiFab(a_ba, a_dm, a_ntracers, a_nghost_force, MFInfo(), a_eb));
    m_tracer_forces->setVal(0.0);
  }

  m_ep_u_mac.reset( new MultiFab(u_ba, a_dm, 1, a_nghost_mac, MFInfo(), a_eb) );
  m_ep_u_mac->setVal(0.0);

  m_ep_v_mac.reset( new MultiFab(v_ba, a_dm, 1, a_nghost_mac, MFInfo(), a_eb) );
  m_ep_v_mac->setVal(0.0);

  m_ep_w_mac.reset( new MultiFab(w_ba, a_dm, 1, a_nghost_mac, MFInfo(), a_eb) );
  m_ep_w_mac->setVal(0.0);

  m_RHS_proj.reset( new MultiFab(a_ba, a_dm, 1, a_nghost_mac, MFInfo(), a_eb) );
  m_RHS_proj->setVal(0.0);

  m_depdt.reset( new MultiFab(a_ba, a_dm, 1, a_nghost_state, MFInfo(), a_eb) );
  m_depdt->setVal(0.0);

  m_rho_nph.reset( new MultiFab(a_ba, a_dm, 1, a_nghost_state, MFInfo(), a_eb) );
  m_rho_nph->setVal(0.0);

  m_epf_nph.reset( new MultiFab(a_ba, a_dm, 1, a_nghost_state, MFInfo(), a_eb) );
  m_epf_nph->setVal(1.0);

  if (a_nspecies > 0) {
    m_J_k[0].reset( new MultiFab(u_ba, a_dm, a_nspecies, 0, MFInfo(), a_eb) );
    m_J_k[0]->setVal(0.0);

    m_J_k[1].reset( new MultiFab(v_ba, a_dm, a_nspecies, 0, MFInfo(), a_eb) );
    m_J_k[1]->setVal(0.0);

    m_J_k[2].reset( new MultiFab(w_ba, a_dm, a_nspecies, 0, MFInfo(), a_eb) );
    m_J_k[2]->setVal(0.0);
  }

  if (a_nspecies > 0 && a_enthalpy) {
    m_h_k[0].reset( new MultiFab(u_ba, a_dm, a_nspecies, 0, MFInfo(), a_eb) );
    m_h_k[0]->setVal(0.0);

    m_h_k[1].reset( new MultiFab(v_ba, a_dm, a_nspecies, 0, MFInfo(), a_eb) );
    m_h_k[1]->setVal(0.0);

    m_h_k[2].reset( new MultiFab(w_ba, a_dm, a_nspecies, 0, MFInfo(), a_eb) );
    m_h_k[2]->setVal(0.0);
  }
}


void StepData::
define ( int const a_nlev,
         int const a_enthalpy,
         int const a_nspecies,
         int const a_ntracers)
{
  m_nlev = a_nlev;

  m_step_data.resize( m_nlev );

  m_has_enthalpy = (a_enthalpy);
  m_has_species  = (a_nspecies >  0) ? 1 : 0;
  m_has_tracer   = (a_ntracers >  0) ? 1 : 0;
}

// StepData Methods that return vectors of constant MultiFabs
//-------------------------------------------------------------------------------------------//

Vector< MultiFab const* > StepData::conv_velocity_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_conv_velocity.get());
  }
  return r;
}


Vector< MultiFab const* > StepData::conv_scalars_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_conv_scalars.get());
  }
  return r;
}

Vector< MultiFab const* > StepData::conv_species_const () const noexcept
{ AMREX_ASSERT( m_nlev > 0 );
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_conv_species.get());
  }
  return r;
}

Vector< MultiFab const* > StepData::lap_T_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_lap_T.get());
  }
  return r;
}

Vector< MultiFab const* > StepData::lap_tracer_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_lap_tracer.get());
  }
  return r;
}

Vector< MultiFab const* > StepData::div_J_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if ( m_has_species ) {
      r.push_back(m_step_data[lev]->m_div_J.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab const* > StepData::div_hJ_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_div_hJ.get());
  }
  return r;
}


// StepData Methods that return vectors of MultiFabs
//-------------------------------------------------------------------------------------------//

Vector< MultiFab    * > StepData::conv_velocity () const noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_conv_velocity.get());
  }
  return r;
}

Vector< MultiFab* > StepData::conv_scalars () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_conv_scalars.get());
  }
  return r;
}

Vector< MultiFab* > StepData::conv_species () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_conv_species.get());
  }
  return r;
}

Vector< MultiFab* > StepData::lap_T () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_lap_T.get());
  }
  return r;
}

Vector< MultiFab* > StepData::lap_tracer () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_lap_tracer.get());
  }
  return r;
}

Vector< MultiFab* > StepData::div_J () const noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if ( m_has_species ) {
      r.push_back(m_step_data[lev]->m_div_J.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab* > StepData::div_hJ () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_step_data[lev]->m_div_hJ.get());
  }
  return r;
}


// StepData Methods that return constant Array4 pointers
//-------------------------------------------------------------------------------------------//

Array4<Real const> StepData::
conv_velocity_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  AMREX_ASSERT( a_scomp < AMREX_SPACEDIM );
  return m_step_data[a_lev]->m_conv_velocity->const_array(a_mfi, a_scomp);
}

Array4<Real const> StepData::
conv_scalars_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  AMREX_ASSERT( a_scomp < m_step_data[a_lev]->m_conv_scalars->nComp() );
  return m_step_data[a_lev]->m_conv_scalars->const_array(a_mfi, a_scomp);
}

Array4<Real const> StepData::
lap_T_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_enthalpy )
      ? m_step_data[a_lev]->m_lap_T->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

Array4<Real const> StepData::
lap_tracer_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_tracer )
      ? m_step_data[a_lev]->m_lap_tracer->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

Array4<Real const> StepData::
conv_species_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept {
  AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_species )
      ? m_step_data[a_lev]->m_conv_species->const_array(a_mfi, a_scomp)
      : Array4<Real const>{} ;
}

Array4<Real const> StepData::
div_J_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept {
  AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_species )
      ? m_step_data[a_lev]->m_div_J->const_array(a_mfi, a_scomp)
      : Array4<Real const>{} ;
}

Array4<Real const> StepData::
div_hJ_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept {
  AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_species && m_has_enthalpy )
      ? m_step_data[a_lev]->m_div_hJ->const_array(a_mfi, a_scomp)
      : Array4<Real const>{} ;
}


// StepData Methods that return Array4 pointers
//-------------------------------------------------------------------------------------------//

Array4<Real      > StepData::
conv_velocity ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  AMREX_ASSERT( a_scomp < AMREX_SPACEDIM );
  return m_step_data[a_lev]->m_conv_velocity->array(a_mfi, a_scomp);
}


Array4<Real      > StepData::
conv_scalars ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  AMREX_ASSERT( a_scomp < m_step_data[a_lev]->m_conv_scalars->nComp() );
  return m_step_data[a_lev]->m_conv_scalars->array(a_mfi, a_scomp);
}

Array4<Real      > StepData::
lap_T ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_enthalpy )
      ? m_step_data[a_lev]->m_lap_T->array(a_mfi, a_scomp)
      : Array4<Real      >{};
}

Array4<Real      > StepData::
lap_tracer ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_tracer )
      ? m_step_data[a_lev]->m_lap_tracer->array(a_mfi, a_scomp)
      : Array4<Real      >{};
}

Array4<Real      > StepData::
conv_species ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_species )
      ? m_step_data[a_lev]->m_conv_species->array(a_mfi, a_scomp)
      : Array4<Real      >{};
}

Array4<Real      > StepData::
div_J ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_species )
      ? m_step_data[a_lev]->m_div_J->array(a_mfi, a_scomp)
      : Array4<Real      >{};
}

Array4<Real      > StepData::
div_hJ ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_nlev );
  return ( m_has_species && m_has_enthalpy )
      ? m_step_data[a_lev]->m_div_hJ->array(a_mfi, a_scomp)
      : Array4<Real      >{};
}




MFIXStepData::
MFIXStepData ( int const a_nlev,
               Vector<BoxArray>                  const& a_ba,
               Vector<DistributionMapping>       const& a_dm,
               Vector<EBFArrayBoxFactory const*> const& a_factory,
               int const a_nghost_state, int const a_nghost_force,
               MFIXBoundaryConditions& a_bcs)
  : m_nlev(a_nlev)
  , m_ba(a_ba)
  , m_dm(a_dm)
  , m_factory(a_factory)
  , m_nghost_state(a_nghost_state)
  , m_nghost_force(a_nghost_force)
  , m_bcs(a_bcs)
{ }


void
MFIXStepData::
define ( Real const a_time,
         int const a_density,
         int const a_enthalpy,
         int const a_nspecies,
         int const a_ntracers,
         int const a_eb_flow,
         int const a_eb_temperature )
{
  m_density  = a_density;
  m_enthalpy = a_enthalpy;
  m_nspecies = a_nspecies;
  m_ntracers = a_ntracers;

  m_ncomp_adv = 2 + m_ntracers;

  m_shared_data.resize( m_nlev );

  m_predictor.define( m_nlev, a_enthalpy, a_nspecies, a_ntracers);
  m_corrector.define( m_nlev, a_enthalpy, a_nspecies, a_ntracers);

  for (int lev(0); lev<m_nlev; ++lev) {

    m_shared_data[lev].reset( new SharedStepData(m_ba[lev], m_dm[lev], *m_factory[lev],
        m_enthalpy, m_nspecies, m_ntracers, m_ncomp_adv,
        m_nghost_state, m_nghost_force, m_nghost_mac) );

    m_predictor.m_step_data[lev].reset( new UniqueStepData(m_ba[lev], m_dm[lev], *m_factory[lev],
        m_enthalpy, m_nspecies, m_ntracers, m_ncomp_adv) );

  }

  if ( a_eb_flow ) { include_eb_flow( a_time ); }
  if ( a_enthalpy && a_eb_temperature ) { include_eb_temperature( a_time ); }

  m_defined = 1;
}


void
MFIXStepData::
include_corrector ()
{ AMREX_ALWAYS_ASSERT( m_defined );

  m_has_corrector = 1;

  for (int lev(0); lev<m_nlev; ++lev) {
    m_corrector.m_step_data[lev].reset( new UniqueStepData(m_ba[lev], m_dm[lev], *m_factory[lev],
        m_enthalpy, m_nspecies, m_ntracers, m_ncomp_adv) );
  }
}


void
MFIXStepData::
include_eb_flow ( Real const a_time )
{
  m_has_eb_flow = 1;

  for (int lev(0); lev<m_nlev; lev++) {

    const BoxArray&              ba = m_ba[lev];
    const DistributionMapping&   dm = m_dm[lev];
    const FabFactory<FArrayBox>& eb = *m_factory[lev];

    m_shared_data[lev]->m_eb_flow_velocity.reset(
        new MultiFab(ba, dm, 3, m_nghost_state, MFInfo(), eb) );

    m_shared_data[lev]->m_eb_flow_velocity->setVal(0.0);

    m_shared_data[lev]->m_eb_flow_scalars.reset(
        new MultiFab(ba, dm, m_ncomp_adv, m_nghost_state, MFInfo(), eb) );

    m_shared_data[lev]->m_eb_flow_scalars->setVal(0.0);

    if (m_nspecies > 0) {
      m_shared_data[lev]->m_eb_flow_species.reset(
          new MultiFab(ba, dm, m_nspecies, m_nghost_state, MFInfo(), eb) );

      m_shared_data[lev]->m_eb_flow_species->setVal(0.0);

    }
  } // nlev

  m_bcs.set_eb_velocity_bcs( a_time, eb_flow_velocity() );

  m_bcs.set_eb_scalar_bcs( eb_flow_scalars(),
      eb_flow_species(), m_enthalpy, m_nspecies );
}


void
MFIXStepData::
include_eb_temperature ( Real const /*a_time*/ )
{
  m_has_eb_temperature = 1;

  for (int lev(0); lev<m_nlev; lev++) {

    const BoxArray&              ba = m_ba[lev];
    const DistributionMapping&   dm = m_dm[lev];
    const FabFactory<FArrayBox>& eb = *m_factory[lev];

    m_shared_data[lev]->m_eb_temperature.reset(
        new MultiFab(ba, dm, 1, m_nghost_state, MFInfo(), eb) );

    m_shared_data[lev]->m_eb_temperature->setVal(0.0);
  }
}


void
MFIXStepData::
include_pressure ()
{ AMREX_ALWAYS_ASSERT( m_defined );

  m_has_pert_pressure = 1;

  for (int lev(0); lev<m_nlev; lev++) {

    const BoxArray&              ba = m_ba[lev];
    const DistributionMapping&   dm = m_dm[lev];
    const FabFactory<FArrayBox>& eb = *m_factory[lev];

    const BoxArray& nd_ba = amrex::convert(ba, IntVect{1,1,1});
    m_shared_data[lev]->m_pert_pressure.reset(
        new MultiFab(nd_ba, dm, 1, m_nghost_state, MFInfo(), eb) );

    m_shared_data[lev]->m_pert_pressure->setVal(0.0);
  }
}


void
MFIXStepData::
set_pressure ( Vector< MultiFab const* > a_pressure )
{ AMREX_ASSERT( m_defined && m_has_pert_pressure );
  AMREX_ASSERT( a_pressure.size() == m_nlev );

  for( int lev(0); lev<m_nlev; ++lev) {

    int const ncomp(a_pressure[lev]->nComp());
    int const ngrow(a_pressure[lev]->nGrow());

    AMREX_ASSERT( ncomp == m_shared_data[lev]->m_pert_pressure->nComp() );
    AMREX_ASSERT( ngrow == m_shared_data[lev]->m_pert_pressure->nGrow() );

    MultiFab::Copy(*m_shared_data[lev]->m_pert_pressure,
        *a_pressure[lev], 0, 0, ncomp, ngrow);
  }
}


void
MFIXStepData::
include_acceleration ()
{ AMREX_ALWAYS_ASSERT( m_defined );

  m_has_acceleration = 1;

  for (int lev(0); lev<m_nlev; lev++) {

    m_shared_data[lev]->m_acceleration.reset( new MultiFab(m_ba[lev], m_dm[lev], 3,
        m_nghost_state, MFInfo(), *m_factory[lev]) );

    m_shared_data[lev]->m_acceleration->setVal(0.0);
  }
}

void
MFIXStepData::
include_avg_particle_data (int ncomp)
{ AMREX_ALWAYS_ASSERT( m_defined );

  if (ncomp > 0) {
    m_has_avg_particle_data = 1;

    for (int lev(0); lev<m_nlev; lev++) {

      m_shared_data[lev]->m_avg_particle_data.reset( new MultiFab(m_ba[lev], m_dm[lev],
          ncomp, 1, MFInfo(), *m_factory[lev]) );

      m_shared_data[lev]->m_avg_particle_data->setVal(0.0);
    }
  }
}



// MFIXStepData Methods that return vectors of constant MultiFabs
//-------------------------------------------------------------------------------------------//

Vector< MultiFab const* > MFIXStepData::vel_S_p_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_vel_S_p.get());
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::vel_S_c_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_vel_S_c.get());
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::vel_forces_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_vel_forces.get());
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::tracer_forces_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_tracer_forces.get());
  }
  return r;
}

Vector<Array<MultiFab const*,3>> MFIXStepData::mac_vel_const () const noexcept {
  Vector<Array<MultiFab const*,3>> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back({AMREX_D_DECL(
        m_shared_data[lev]->m_ep_u_mac.get(),
        m_shared_data[lev]->m_ep_v_mac.get(),
        m_shared_data[lev]->m_ep_w_mac.get() )});
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::ep_u_mac_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_ep_u_mac.get());
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::ep_v_mac_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_ep_v_mac.get());
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::ep_w_mac_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_ep_w_mac.get());
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::RHS_proj_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_RHS_proj.get());
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::depdt_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_depdt.get());
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::epf_nph_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_epf_nph.get());
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::rho_nph_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_rho_nph.get());
  }
  return r;
}

Vector< Array< MultiFab const*, AMREX_SPACEDIM>> MFIXStepData::J_k_const () const noexcept
{
  Vector<Array<MultiFab const*, AMREX_SPACEDIM>> r;
  r.reserve(m_nlev);
  if (m_nspecies > 0 ) {
    for (int lev(0); lev < m_nlev; ++lev) {
      r.push_back({AMREX_D_DECL(
          m_shared_data[lev]->m_J_k[0].get(),
          m_shared_data[lev]->m_J_k[1].get(),
          m_shared_data[lev]->m_J_k[2].get() )});
    }
  }
  return r;
}

Vector< Array< MultiFab const*, AMREX_SPACEDIM>> MFIXStepData::h_k_const () const noexcept
{
  Vector<Array<MultiFab const*, AMREX_SPACEDIM>> r;
  r.reserve(m_nlev);
  if (m_nspecies > 0 && m_enthalpy) {
    for (int lev(0); lev < m_nlev; ++lev) {
      r.push_back({AMREX_D_DECL(
          m_shared_data[lev]->m_h_k[0].get(),
          m_shared_data[lev]->m_h_k[1].get(),
          m_shared_data[lev]->m_h_k[2].get() )});
    }
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::eb_flow_velocity_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_eb_flow ) {
      r.push_back(m_shared_data[lev]->m_eb_flow_velocity.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::eb_flow_scalars_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_eb_flow ) {
      r.push_back(m_shared_data[lev]->m_eb_flow_scalars.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::eb_flow_species_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_eb_flow && (m_nspecies > 0)) {
      r.push_back(m_shared_data[lev]->m_eb_flow_species.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::eb_temperature_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_eb_temperature) {
      r.push_back(m_shared_data[lev]->m_eb_temperature.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::pert_pressure_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_pert_pressure) {
      r.push_back(m_shared_data[lev]->m_pert_pressure.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::acceleration_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_acceleration) {
      r.push_back(m_shared_data[lev]->m_acceleration.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab const* > MFIXStepData::avg_particle_data_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_avg_particle_data) {
      r.push_back(m_shared_data[lev]->m_avg_particle_data.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}


// MFIXStepData Methods that return vectors of MultiFabs
//-------------------------------------------------------------------------------------------//

Vector<Array<MultiFab*,3>> MFIXStepData::mac_vel () const noexcept {
  Vector<Array<MultiFab*,3>> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back({AMREX_D_DECL(
        m_shared_data[lev]->m_ep_u_mac.get(),
        m_shared_data[lev]->m_ep_v_mac.get(),
        m_shared_data[lev]->m_ep_w_mac.get() )});
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::ep_u_mac () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_ep_u_mac.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::ep_v_mac () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_ep_v_mac.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::ep_w_mac () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_ep_w_mac.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::RHS_proj () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_RHS_proj.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::depdt () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_depdt.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData:: epf_nph () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_epf_nph.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::rho_nph () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_rho_nph.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::vel_S_p () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_vel_S_p.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::vel_S_c () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_vel_S_c.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::vel_forces () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_vel_forces.get());
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::tracer_forces () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    r.push_back(m_shared_data[lev]->m_tracer_forces.get());
  }
  return r;
}

Vector< Array< MultiFab*, AMREX_SPACEDIM>> MFIXStepData::J_k () const noexcept
{
  Vector<Array<MultiFab*, AMREX_SPACEDIM>> r;
  r.reserve(m_nlev);
  if (m_nspecies > 0 ) {
    for (int lev(0); lev < m_nlev; ++lev) {
      r.push_back({AMREX_D_DECL(
          m_shared_data[lev]->m_J_k[0].get(),
          m_shared_data[lev]->m_J_k[1].get(),
          m_shared_data[lev]->m_J_k[2].get() )});
    }
  }
  return r;
}

Vector< Array< MultiFab*, AMREX_SPACEDIM>> MFIXStepData::h_k () const noexcept
{
  Vector<Array<MultiFab*, AMREX_SPACEDIM>> r;
  r.reserve(m_nlev);
  if (m_nspecies > 0 && m_enthalpy) {
    for (int lev(0); lev < m_nlev; ++lev) {
      r.push_back({AMREX_D_DECL(
          m_shared_data[lev]->m_h_k[0].get(),
          m_shared_data[lev]->m_h_k[1].get(),
          m_shared_data[lev]->m_h_k[2].get() )});
    }
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::eb_flow_velocity () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_eb_flow ) {
      r.push_back(m_shared_data[lev]->m_eb_flow_velocity.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::eb_flow_scalars () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_eb_flow ) {
      r.push_back(m_shared_data[lev]->m_eb_flow_scalars.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::eb_flow_species () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_eb_flow && (m_nspecies > 0)) {
      r.push_back(m_shared_data[lev]->m_eb_flow_species.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::eb_temperature () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_eb_temperature) {
      r.push_back(m_shared_data[lev]->m_eb_temperature.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::pert_pressure () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_pert_pressure) {
      r.push_back(m_shared_data[lev]->m_pert_pressure.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::acceleration () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_acceleration) {
      r.push_back(m_shared_data[lev]->m_acceleration.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}

Vector< MultiFab* > MFIXStepData::avg_particle_data () const noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_nlev);
  for (int lev(0); lev < m_nlev; ++lev) {
    if (m_has_avg_particle_data) {
      r.push_back(m_shared_data[lev]->m_avg_particle_data.get());
    } else {
      r.push_back(nullptr);
    }
  }
  return r;
}


// MFIXStepData Methods that return constant MultiFab pointers
//-------------------------------------------------------------------------------------------//

Array<MultiFab const*,3> MFIXStepData::mac_vel_const ( int const a_lev ) const noexcept {
  return { m_shared_data[a_lev]->m_ep_u_mac.get(),
           m_shared_data[a_lev]->m_ep_v_mac.get(),
           m_shared_data[a_lev]->m_ep_w_mac.get() };
}

MultiFab const* MFIXStepData::depdt_const ( int const a_lev ) const noexcept
{ return m_shared_data[a_lev]->m_depdt.get(); }

MultiFab const* MFIXStepData::vel_S_p_const ( int const a_lev ) const noexcept
{ return m_shared_data[a_lev]->m_vel_S_p.get(); }

MultiFab const* MFIXStepData::vel_S_c_const ( int const a_lev ) const noexcept
{ return m_shared_data[a_lev]->m_vel_S_c.get(); }


// MFIXStepData Methods that return MultiFab pointers
//-------------------------------------------------------------------------------------------//

Array<MultiFab      *,3> MFIXStepData::mac_vel ( int const a_lev ) const noexcept {
  return { m_shared_data[a_lev]->m_ep_u_mac.get(),
           m_shared_data[a_lev]->m_ep_v_mac.get(),
           m_shared_data[a_lev]->m_ep_w_mac.get() };
}

MultiFab      * MFIXStepData::depdt ( int const a_lev ) const noexcept
{ return m_shared_data[a_lev]->m_depdt.get(); }

MultiFab      * MFIXStepData::epf_nph ( int const a_lev ) const noexcept
{ return m_shared_data[a_lev]->m_epf_nph.get(); }

MultiFab      * MFIXStepData::rho_nph ( int const a_lev ) const noexcept
{ return m_shared_data[a_lev]->m_rho_nph.get(); }


// MFIXStepData Methods that return constant Array4 pointers
//-------------------------------------------------------------------------------------------//

Array4<Real const> MFIXStepData::
depdt_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_depdt->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXStepData::
epf_nph_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_epf_nph->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXStepData::
rho_nph_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_rho_nph->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXStepData::
vel_S_p_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_vel_S_p->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXStepData::
vel_S_c_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_vel_S_c->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXStepData::
vel_forces_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_vel_forces->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXStepData::
tracer_forces_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_tracer_forces->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXStepData::
acceleration_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_acceleration->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXStepData::
avg_particle_data_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_avg_particle_data->const_array(a_mfi, a_scomp);
}


// MFIXStepData Methods that return Array4 pointers
//-------------------------------------------------------------------------------------------//

Array4<Real      > MFIXStepData::
depdt ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_depdt->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXStepData::
epf_nph ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_epf_nph->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXStepData::
rho_nph ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_rho_nph->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXStepData::
vel_S_p ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_vel_S_p->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXStepData::
vel_S_c ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_vel_S_c->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXStepData::
vel_forces ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_vel_forces->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXStepData::
tracer_forces ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_tracer_forces->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXStepData::
acceleration ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_acceleration->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXStepData::
avg_particle_data ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{ AMREX_ASSERT( a_lev < m_shared_data.size() );
  return m_shared_data[a_lev]->m_avg_particle_data->array(a_mfi, a_scomp);
}
