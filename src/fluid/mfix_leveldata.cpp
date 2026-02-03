#include <AMReX_EBMultiFabUtil.H>

#include <mfix_leveldata.H>

using namespace amrex;

void LevelData::
define ( int const a_nghost,
         BoxArray const& a_ba,
         DistributionMapping const& a_dmap,
         FabFactory<FArrayBox> const& a_factory)
{
  AMREX_ASSERT( !defined() );

  // Fluid volume fraction (void fraction)
  m_epf_old.reset( new MultiFab(a_ba, a_dmap, 1, a_nghost, MFInfo(), a_factory) );
  m_epf.reset( new MultiFab(a_ba, a_dmap, 1, a_nghost, MFInfo(), a_factory) );

  // Fluid density
  m_rho_old.reset( new MultiFab(a_ba, a_dmap, 1, a_nghost, MFInfo(), a_factory) );
  m_rho.reset( new MultiFab(a_ba, a_dmap, 1, a_nghost, MFInfo(), a_factory) );

  // Fluid velocity
  m_vel_old.reset( new MultiFab(a_ba, a_dmap, 3, a_nghost, MFInfo(), a_factory) );
  m_vel.reset( new MultiFab(a_ba, a_dmap, 3, a_nghost, MFInfo(), a_factory) );

  // perturbational pressure
  m_pert_p.reset( new MultiFab(amrex::convert(a_ba, IntVect{1,1,1}), a_dmap,
      /*ncomp =*/1, a_nghost, MFInfo(), a_factory) );

  // gradient of perturbational pressure
  m_grad_p.reset( new MultiFab(a_ba, a_dmap, /*ncomp =*/3, a_nghost, MFInfo(), a_factory));

  // Fluid temperature
  if ( has_temperature() ) {
    m_T_old.reset( new MultiFab(a_ba, a_dmap, 1, a_nghost, MFInfo(), a_factory) );
    m_T.reset( new MultiFab(a_ba, a_dmap, 1, a_nghost, MFInfo(), a_factory) );
  }

  // Fluid (mixture) enthalpy
  if ( has_enthalpy() ) {
    m_h_old.reset( new MultiFab(a_ba, a_dmap, 1, a_nghost, MFInfo(), a_factory) );
    m_h.reset( new MultiFab(a_ba, a_dmap, 1, a_nghost, MFInfo(), a_factory) );
  }

  // Fluid species mass fractions
  if ( has_species() ) {
    int const ncomp( m_fluid.nspecies() );
    m_X_old.reset( new MultiFab(a_ba, a_dmap, ncomp, a_nghost, MFInfo(), a_factory) );
    m_X.reset( new MultiFab(a_ba, a_dmap, ncomp, a_nghost, MFInfo(), a_factory) );
  }

  // Fluid tracers (passive)
  if ( has_tracer() ) {
    int const ncomp( m_fluid.ntracer() );
    m_tracer_old.reset( new MultiFab(a_ba, a_dmap, ncomp, a_nghost, MFInfo(), a_factory) );
    m_tracer.reset( new MultiFab(a_ba, a_dmap, ncomp, a_nghost, MFInfo(), a_factory) );
  }

  // Interphase transfer (drag, heat, etc)
  { InterphaseTxfrIndexes txfr_comps;
    int const ncomp(txfr_comps.count);
    m_txfr.reset( new MultiFab(a_ba, a_dmap, ncomp, a_nghost, MFInfo(), a_factory) );
  }

  // Interphase chemical reaction source terms
  { InterphaseChemTxfrIndexes chem_txfr_comps(m_fluid.nspecies(), m_rxns.solve());
    int const ncomp(chem_txfr_comps.count);
    m_chem_txfr.reset( new MultiFab(a_ba, a_dmap, ncomp, a_nghost, MFInfo(), a_factory));
  }

  m_vorticity.reset( new MultiFab(a_ba, a_dmap, /*ncomp =*/1, a_nghost, MFInfo(), a_factory) );
  m_mac_phi.reset( new MultiFab(a_ba, a_dmap, /*ncomp =*/1, a_nghost, MFInfo(), a_factory) );
  m_div_tau.reset( new MultiFab(a_ba, a_dmap, /*ncomp =*/3, 0, MFInfo(), a_factory) );

  m_defined = 1;
}


void LevelData::
clear (  )
{
  if ( defined() ) {
    // Fluid volume fraction (void fraction)
    m_epf_old.reset( nullptr );
    m_epf.reset( nullptr );

    // Fluid density
    m_rho_old.reset( nullptr );
    m_rho.reset( nullptr );

    // Fluid velocity
    m_vel_old.reset( nullptr );
    m_vel.reset( nullptr );

    // perturbational pressure
    m_pert_p.reset( nullptr );

    // gradient of perturbational pressure
    m_grad_p.reset( nullptr );

    // Fluid temperature
    if ( has_temperature() ) {
      m_T_old.reset( nullptr );
      m_T.reset( nullptr );
    }

    // Fluid (mixture) enthalpy
    if ( has_enthalpy() ) {
      m_h_old.reset( nullptr );
      m_h.reset( nullptr );
    }

    // Fluid species mass fractions
    if ( has_species() ) {
      m_X_old.reset( nullptr );
      m_X.reset( nullptr );
    }

    // Fluid tracers (passive)
    if ( has_tracer() ) {
      m_tracer_old.reset( nullptr );
      m_tracer.reset( nullptr );
    }

    // Interphase transfer (drag, heat, etc)
    m_txfr.reset( nullptr );

    // Interphase chemical reaction source terms
    m_chem_txfr.reset( nullptr );

    m_vorticity.reset( nullptr );
    m_mac_phi.reset( nullptr );
    m_div_tau.reset( nullptr );
  }
  m_defined = 0;
}

void LevelData::
initVals ( Real const a_value)
{ AMREX_ASSERT( defined() );

  // Fluid volume fraction (void fraction)
  m_epf_old->setVal( 1.0 );
  m_epf->setVal( 1.0 );

  // Fluid density
  m_rho_old->setVal( a_value );
  m_rho->setVal( a_value );

  // Fluid velocity
  m_vel_old->setVal( a_value );
  m_vel->setVal( a_value );

  // perturbational pressure and gradient
  m_pert_p->setVal( a_value );
  m_grad_p->setVal( 0. );

  // Fluid temperature
  if ( has_temperature() ) {
    m_T_old->setVal( a_value );
    m_T->setVal( a_value );
  }

  // Fluid (mixture) enthalpy
  if ( has_enthalpy() ) {
    m_h_old->setVal( a_value );
    m_h->setVal( a_value );
  }

  // Fluid species mass fractions
  if ( has_species() ) {
    m_X_old->setVal( a_value );
    m_X->setVal( a_value );
  }

  // Fluid tracers (passive)
  if ( has_tracer() ) {
    m_tracer_old->setVal( a_value );
    m_tracer->setVal( a_value );
  }

  // Interphase transfer (drag, heat, etc) and chemistry sources
  m_txfr->setVal( 0. );
  m_chem_txfr->setVal( 0. );

  m_vorticity->setVal( a_value );
  m_mac_phi->setVal( 0. );
  m_div_tau->setVal( a_value );
}



void LevelData::
swapOldAndNew () {
  AMREX_ASSERT( defined() );

  // Fluid volume fraction (void fraction)
  //m_epf_old.swap( m_epf );

  // Fluid density
  m_rho_old.swap( m_rho );

  // Fluid velocity
  m_vel_old.swap( m_vel );

  // Fluid temperature
  if ( has_temperature() ) {
    m_T_old.swap( m_T );
  }

  // Fluid (mixture) enthalpy
  if ( has_enthalpy() ) {
    m_h_old.swap( m_h );
  }

  // Fluid species mass fractions
  if ( has_species() ) {
    m_X_old.swap( m_X );
  }

  // Fluid tracers (passive)
  if ( has_tracer() ) {
    m_tracer_old.swap( m_tracer );
  }
}


void LevelData::
resetNewWithOld () {
  AMREX_ASSERT( defined() );

  MultiFab::Copy(*m_vel, *m_vel_old, 0, 0, m_vel->nComp(), m_vel->nGrow());

  if (m_fluid.solve_density()) {
    MultiFab::Copy(*m_rho, *m_rho_old, 0, 0, m_rho->nComp(), m_rho->nGrow());
  }

  if ( has_temperature() ) {
    MultiFab::Copy(*m_T, *m_T_old, 0, 0, m_T->nComp(), m_T->nGrow());
  }

  if ( has_enthalpy() ) {
    MultiFab::Copy(*m_h, *m_h_old, 0, 0, m_h->nComp(), m_h->nGrow());
  }

  if ( has_tracer() ) {
    MultiFab::Copy(*m_tracer, *m_tracer_old, 0, 0, m_tracer->nComp(), m_tracer->nGrow());
  }

  if ( has_species() ) {
    MultiFab::Copy(*m_X, *m_X_old, 0, 0, m_X->nComp(), m_X->nGrow());
  }
}


void LevelData::
resetOldWithNew () {
  AMREX_ASSERT( defined() );

  MultiFab::Copy(*m_vel_old, *m_vel, 0, 0, m_vel->nComp(), m_vel->nGrow());

  MultiFab::Copy(*m_rho_old, *m_rho, 0, 0, m_rho->nComp(), m_rho->nGrow());

  if ( has_temperature() ) {
    MultiFab::Copy(*m_T_old, *m_T, 0, 0, m_T->nComp(), m_T->nGrow());
  }

  if ( has_enthalpy() ) {
    MultiFab::Copy(*m_h_old, *m_h, 0, 0, m_h->nComp(), m_h->nGrow());
  }

  if ( has_tracer() ) {
    MultiFab::Copy(*m_tracer_old, *m_tracer, 0, 0, m_tracer->nComp(), m_tracer->nGrow());
  }

  if ( has_species() ) {
    MultiFab::Copy(*m_X_old, *m_X, 0, 0, m_X->nComp(), m_X->nGrow());
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////



//  OLD  ////////////////////////////////////////////////////////////////////////////////////////////

Vector<MultiFab      *> MFIXLevelData::epf_old () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_epf_old.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::rho_old () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    r.push_back(m_level_data[lev]->m_rho_old.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::vel_old () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_vel_old.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::T_old () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_T_old.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::h_old () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_h_old.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::X_old () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_X_old.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::tracer_old () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_tracer_old.get());
  }
  return r;
}


//  NEW  ////////////////////////////////////////////////////////////////////////////////////////////

Vector<MultiFab      *> MFIXLevelData::epf () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_epf.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::rho () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_rho.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::vel () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_vel.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::T () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_T.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::h () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_h.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::X () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_X.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::tracer () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_tracer.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::pert_p () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_pert_p.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::grad_p () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_grad_p.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::txfr () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_txfr.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::chem_txfr () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_chem_txfr.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::vorticity () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_vorticity.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::div_tau () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_div_tau.get());
  }
  return r;
}

Vector<MultiFab      *> MFIXLevelData::mac_phi () noexcept
{
  Vector<MultiFab      *> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_mac_phi.get());
  }
  return r;
}




//  OLD CONST  //////////////////////////////////////////////////////////////////////////////////////

Vector<MultiFab const*> MFIXLevelData::epf_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_epf_old.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::rho_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_rho_old.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::vel_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_vel_old.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::T_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_T_old.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::h_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_h_old.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::X_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_X_old.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::tracer_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_tracer_old.get());
  }
  return r;
}


//  NEW CONST  //////////////////////////////////////////////////////////////////////////////////////

Vector<MultiFab const*> MFIXLevelData::epf_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_epf.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::rho_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_rho.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::vel_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_vel.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::T_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_T.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::h_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_h.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::X_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_X.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::tracer_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_tracer.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::pert_p_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_pert_p.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::grad_p_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_grad_p.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::txfr_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_txfr.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::chem_txfr_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_chem_txfr.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::vorticity_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_vorticity.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::div_tau_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_div_tau.get());
  }
  return r;
}

Vector<MultiFab const*> MFIXLevelData::mac_phi_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_finest_level+1);
  for (int lev(0); lev <= m_finest_level; ++lev) {
    AMREX_ASSERT( m_level_data[lev]->defined() );
    r.push_back(m_level_data[lev]->m_mac_phi.get());
  }
  return r;
}


//  OLD LEV /////////////////////////////////////////////////////////////////////////////////////////

MultiFab      * MFIXLevelData::epf_old ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_epf_old.get();
}

MultiFab      * MFIXLevelData::rho_old ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_rho_old.get();
}

MultiFab      * MFIXLevelData::vel_old ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vel_old.get();
}

MultiFab      * MFIXLevelData::T_old ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_T_old.get();
}

MultiFab      * MFIXLevelData::h_old ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_h_old.get();
}

MultiFab      * MFIXLevelData::X_old ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_X_old.get();
}

MultiFab      * MFIXLevelData::tracer_old ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_tracer_old.get();
}


//  NEW LEV /////////////////////////////////////////////////////////////////////////////////////////

MultiFab      * MFIXLevelData::epf ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_epf.get();
}

MultiFab      * MFIXLevelData::rho ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_rho.get();
}

MultiFab      * MFIXLevelData::vel ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vel.get();
}

MultiFab      * MFIXLevelData::T ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_T.get();
}

MultiFab      * MFIXLevelData::h ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_h.get();
}

MultiFab      * MFIXLevelData::X ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_X.get();
}

MultiFab      * MFIXLevelData::tracer ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_tracer.get();
}

MultiFab      * MFIXLevelData::pert_p ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_pert_p.get();
}

MultiFab      * MFIXLevelData::grad_p ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_grad_p.get();
}

MultiFab      * MFIXLevelData::txfr ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_txfr.get();
}

MultiFab      * MFIXLevelData::chem_txfr ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_chem_txfr.get();
}

MultiFab      * MFIXLevelData::vorticity ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vorticity.get();
}

MultiFab      * MFIXLevelData::div_tau ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_div_tau.get();
}

MultiFab      * MFIXLevelData::mac_phi ( int const a_lev ) noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_mac_phi.get();
}


//  OLD CONST LEV ///////////////////////////////////////////////////////////////////////////////////

MultiFab const* MFIXLevelData::epf_old_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_epf_old.get();
}

MultiFab const* MFIXLevelData::rho_old_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_rho_old.get();
}

MultiFab const* MFIXLevelData::vel_old_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vel_old.get();
}

MultiFab const* MFIXLevelData::T_old_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_T_old.get();
}

MultiFab const* MFIXLevelData::h_old_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_h_old.get();
}

MultiFab const* MFIXLevelData::X_old_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_X_old.get();
}

MultiFab const* MFIXLevelData::tracer_old_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_tracer_old.get();
}


//  NEW CONST LEV ///////////////////////////////////////////////////////////////////////////////////

MultiFab const* MFIXLevelData::epf_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_epf.get();
}

MultiFab const* MFIXLevelData::rho_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_rho.get();
}

MultiFab const* MFIXLevelData::vel_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vel.get();
}

MultiFab const* MFIXLevelData::T_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_T.get();
}

MultiFab const* MFIXLevelData::h_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_h.get();
}

MultiFab const* MFIXLevelData::X_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_X.get();
}

MultiFab const* MFIXLevelData::tracer_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_tracer.get();
}

MultiFab const* MFIXLevelData::pert_p_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_pert_p.get();
}

MultiFab const* MFIXLevelData::grad_p_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_grad_p.get();
}

MultiFab const* MFIXLevelData::txfr_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_txfr.get();
}

MultiFab const* MFIXLevelData::chem_txfr_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_chem_txfr.get();
}

MultiFab const* MFIXLevelData::vorticity_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vorticity.get();
}

MultiFab const* MFIXLevelData::div_tau_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_div_tau.get();
}

MultiFab const* MFIXLevelData::mac_phi_const ( int const a_lev ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_mac_phi.get();
}

//  OLD ARRAY4 //////////////////////////////////////////////////////////////////////////////////////

Array4<Real      > MFIXLevelData::
epf_old ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_epf_old->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
rho_old ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_rho_old->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
vel_old ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vel_old->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
T_old ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_temperature() )
      ? m_level_data[a_lev]->m_T_old->array(a_mfi, a_scomp)
      : Array4<Real      >{};
}

Array4<Real      > MFIXLevelData::
h_old ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_enthalpy() )
      ? m_level_data[a_lev]->m_h_old->array(a_mfi, a_scomp)
      : Array4<Real      >{};
}

Array4<Real      > MFIXLevelData::
X_old ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_species() )
      ? m_level_data[a_lev]->m_X_old->array(a_mfi, a_scomp)
      : Array4<Real      >{};
}

Array4<Real      > MFIXLevelData::
tracer_old ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_tracer() )
      ? m_level_data[a_lev]->m_tracer_old->array(a_mfi, a_scomp)
      : Array4<Real      >{};
}

//  NEW ARRAY4 //////////////////////////////////////////////////////////////////////////////////////

Array4<Real      > MFIXLevelData::
epf ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_epf->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
rho ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_rho->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
vel ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vel->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
T ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_temperature() )
      ? m_level_data[a_lev]->m_T->array(a_mfi, a_scomp)
      : Array4<Real      >{} ;
}

Array4<Real      > MFIXLevelData::
h ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_enthalpy() )
      ? m_level_data[a_lev]->m_h->array(a_mfi, a_scomp)
      : Array4<Real      >{} ;
}

Array4<Real      > MFIXLevelData::
X ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_species() )
      ? m_level_data[a_lev]->m_X->array(a_mfi, a_scomp)
      : Array4<Real      >{} ;
}

Array4<Real      > MFIXLevelData::
tracer ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_tracer() )
      ? m_level_data[a_lev]->m_tracer->array(a_mfi, a_scomp)
      : Array4<Real      >{} ;
}

Array4<Real      > MFIXLevelData::
pert_p ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_pert_p->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
grad_p ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_grad_p->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
txfr ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_txfr->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
chem_txfr ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_chem_txfr->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
vorticity ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vorticity->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
div_tau ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_div_tau->array(a_mfi, a_scomp);
}

Array4<Real      > MFIXLevelData::
mac_phi ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_mac_phi->array(a_mfi, a_scomp);
}

//  OLD CONST ARRAY4 ////////////////////////////////////////////////////////////////////////////////

Array4<Real const> MFIXLevelData::
epf_old_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_epf_old->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
rho_old_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_rho_old->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
vel_old_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vel_old->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
T_old_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_temperature() )
      ? m_level_data[a_lev]->m_T_old->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

Array4<Real const> MFIXLevelData::
h_old_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_enthalpy() )
      ? m_level_data[a_lev]->m_h_old->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

Array4<Real const> MFIXLevelData::
X_old_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_species() )
      ? m_level_data[a_lev]->m_X_old->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

Array4<Real const> MFIXLevelData::
tracer_old_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_tracer() )
      ? m_level_data[a_lev]->m_tracer_old->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

//  NEW ARRAY4 //////////////////////////////////////////////////////////////////////////////////////

Array4<Real const> MFIXLevelData::
epf_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_epf->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
rho_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_rho->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
vel_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vel->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
T_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_temperature() )
      ? m_level_data[a_lev]->m_T->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

Array4<Real const> MFIXLevelData::
h_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_enthalpy() )
      ? m_level_data[a_lev]->m_h->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

Array4<Real const> MFIXLevelData::
X_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_species() )
      ? m_level_data[a_lev]->m_X->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

Array4<Real const> MFIXLevelData::
tracer_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return ( m_level_data[a_lev]->has_tracer() )
      ? m_level_data[a_lev]->m_tracer->const_array(a_mfi, a_scomp)
      : Array4<Real const>{};
}

Array4<Real const> MFIXLevelData::
pert_p_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_pert_p->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
grad_p_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_grad_p->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
txfr_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_txfr->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
chem_txfr_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_chem_txfr->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
vorticity_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_vorticity->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
div_tau_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_div_tau->const_array(a_mfi, a_scomp);
}

Array4<Real const> MFIXLevelData::
mac_phi_const ( int const a_lev, const MFIter& a_mfi, int const a_scomp ) const noexcept
{
  AMREX_ASSERT( m_level_data[a_lev]->defined() );
  return m_level_data[a_lev]->m_mac_phi->const_array(a_mfi, a_scomp);
}
