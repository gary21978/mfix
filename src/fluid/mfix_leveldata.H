#ifndef _MFIX_LEVELDATA_H_
#define _MFIX_LEVELDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MultiFab.H>

#include <mfix_bc.H>
#include <mfix_eb_parms.H>
#include <mfix_fluid.H>
#include <mfix_reactions.H>


class LevelData {

  public:
    // Default Constructor
    LevelData ( MFIXFluidPhase const& a_fluid,
                MFIXReactions  const& a_rxns)
    : m_fluid(a_fluid) , m_rxns(a_rxns)
    { }

    // Destructor
    ~LevelData () = default;

    void define ( int const a_nghost,
                  amrex::BoxArray const& a_ba,
                  amrex::DistributionMapping const& a_dm,
                  amrex::FabFactory<amrex::FArrayBox> const& a_factory);

    void clear ( );

    // Reset initial values
    void initVals ( amrex::Real const a_value = std::numeric_limits<amrex::Real>::max() );

    // Copy old values into new
    void resetNewWithOld ();

    // Copy new values into old
    void resetOldWithNew ();

    // swapp the old and new pointers
    void swapOldAndNew ();

    bool defined () const noexcept { return (m_defined == 1); }

    // Fluid volume fraction (void fraction)
    std::unique_ptr<amrex::MultiFab> m_epf_old;
    std::unique_ptr<amrex::MultiFab> m_epf;

    // Fluid density
    std::unique_ptr<amrex::MultiFab> m_rho_old;
    std::unique_ptr<amrex::MultiFab> m_rho;

    // Fluid velocity
    std::unique_ptr<amrex::MultiFab> m_vel_old;
    std::unique_ptr<amrex::MultiFab> m_vel;

    // Fluid temperature
    std::unique_ptr<amrex::MultiFab> m_T_old;
    std::unique_ptr<amrex::MultiFab> m_T;

    // Fluid (mixture) enthalpy
    std::unique_ptr<amrex::MultiFab> m_h_old;
    std::unique_ptr<amrex::MultiFab> m_h;

    // Fluid species mass fractions
    std::unique_ptr<amrex::MultiFab> m_X_old;
    std::unique_ptr<amrex::MultiFab> m_X;

    // Fluid tracers (passive)
    std::unique_ptr<amrex::MultiFab> m_tracer_old;
    std::unique_ptr<amrex::MultiFab> m_tracer;

    // Fluid perturbational pressure and gradient
    std::unique_ptr<amrex::MultiFab> m_pert_p;
    std::unique_ptr<amrex::MultiFab> m_grad_p;

    // Interphase transfer (drag, heat, etc)
    std::unique_ptr<amrex::MultiFab> m_txfr;
    std::unique_ptr<amrex::MultiFab> m_chem_txfr;

    // Cell-based
    std::unique_ptr<amrex::MultiFab> m_vorticity;
    std::unique_ptr<amrex::MultiFab> m_div_tau;
    std::unique_ptr<amrex::MultiFab> m_mac_phi;

    bool has_temperature () const noexcept
    { return (m_fluid.solve_enthalpy() || m_fluid.constraint.isIdealGas()); }

    bool has_enthalpy () const noexcept
    { return m_fluid.solve_enthalpy(); }

    bool has_species () const noexcept
    { return m_fluid.solve_species(); }

    bool has_tracer () const noexcept
    { return m_fluid.solve_tracer(); }

    MFIXFluidPhase const& fluid () const noexcept { return m_fluid; }

    MFIXReactions const& reactions () const noexcept { return m_rxns; }

  private:

     MFIXFluidPhase const& m_fluid;
     MFIXReactions const&  m_rxns;

     int m_defined = 0;
};

class MFIXLevelData {

  public:

    MFIXLevelData ( int const a_max_level )
    { m_level_data.resize( a_max_level+1 ); }

    ~MFIXLevelData () = default;

    int  finestLevel () const { return m_finest_level; }
    int& finestLevel () { return m_finest_level; }

    void define ( int const a_lev, int const a_nghost,
                  amrex::BoxArray const& a_ba,
                  amrex::DistributionMapping const& a_dm,
                  amrex::FabFactory<amrex::FArrayBox> const& a_factory)
    {
      m_finest_level++;

      AMREX_ASSERT( a_lev == m_finest_level );
      AMREX_ASSERT( !m_level_data[a_lev]->defined() );

      m_level_data[a_lev]->define(a_nghost, a_ba, a_dm, a_factory);
    }

    void move ( int const a_lev, std::unique_ptr<LevelData> a_new_level )
    {
      // This is a new level that we're adding.
      if (a_lev > m_finest_level) {
        AMREX_ASSERT( !m_level_data[a_lev]->defined() );
        m_finest_level++;

      // This is an existing level that we're replacing.
      } else {
        AMREX_ASSERT( m_level_data[a_lev]->defined() );
      }
      m_level_data[a_lev] = std::move( a_new_level);
    }

    void clear ( int const a_lev )
    {
      AMREX_ASSERT( m_level_data[a_lev]->defined() );
      AMREX_ASSERT( a_lev == m_finest_level );

      m_level_data[a_lev].reset();
      m_finest_level--;
    }

    amrex::Vector<amrex::MultiFab      *> epf_old () noexcept;
    amrex::Vector<amrex::MultiFab      *> rho_old () noexcept;
    amrex::Vector<amrex::MultiFab      *> vel_old () noexcept;
    amrex::Vector<amrex::MultiFab      *> T_old () noexcept;
    amrex::Vector<amrex::MultiFab      *> h_old () noexcept;
    amrex::Vector<amrex::MultiFab      *> X_old () noexcept;
    amrex::Vector<amrex::MultiFab      *> tracer_old () noexcept;

    amrex::Vector<amrex::MultiFab      *> epf () noexcept;
    amrex::Vector<amrex::MultiFab      *> rho () noexcept;
    amrex::Vector<amrex::MultiFab      *> vel () noexcept;
    amrex::Vector<amrex::MultiFab      *> T () noexcept;
    amrex::Vector<amrex::MultiFab      *> h () noexcept;
    amrex::Vector<amrex::MultiFab      *> X () noexcept;
    amrex::Vector<amrex::MultiFab      *> tracer () noexcept;

    amrex::Vector<amrex::MultiFab      *> pert_p () noexcept;
    amrex::Vector<amrex::MultiFab      *> grad_p () noexcept;

    amrex::Vector<amrex::MultiFab      *> txfr () noexcept;
    amrex::Vector<amrex::MultiFab      *> chem_txfr () noexcept;

    amrex::Vector<amrex::MultiFab      *> vorticity () noexcept;
    amrex::Vector<amrex::MultiFab      *> div_tau () noexcept;
    amrex::Vector<amrex::MultiFab      *> mac_phi () noexcept;

    amrex::Vector<amrex::MultiFab const*> epf_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> rho_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> vel_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> T_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> h_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> X_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> tracer_old_const () const noexcept;

    amrex::Vector<amrex::MultiFab const*> epf_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> rho_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> vel_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> T_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> h_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> X_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> tracer_const () const noexcept;

    amrex::Vector<amrex::MultiFab const*> pert_p_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> grad_p_const () const noexcept;

    amrex::Vector<amrex::MultiFab const*> txfr_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> chem_txfr_const () const noexcept;

    amrex::Vector<amrex::MultiFab const*> vorticity_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> div_tau_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> mac_phi_const () const noexcept;

    amrex::MultiFab      * epf_old ( int const a_lev ) noexcept;
    amrex::MultiFab      * rho_old ( int const a_lev ) noexcept;
    amrex::MultiFab      * vel_old ( int const a_lev ) noexcept;
    amrex::MultiFab      * T_old ( int const a_lev ) noexcept;
    amrex::MultiFab      * h_old ( int const a_lev ) noexcept;
    amrex::MultiFab      * X_old ( int const a_lev ) noexcept;
    amrex::MultiFab      * tracer_old ( int const a_lev ) noexcept;

    amrex::MultiFab      * epf ( int const a_lev ) noexcept;
    amrex::MultiFab      * rho ( int const a_lev ) noexcept;
    amrex::MultiFab      * vel ( int const a_lev ) noexcept;
    amrex::MultiFab      * T ( int const a_lev ) noexcept;
    amrex::MultiFab      * h ( int const a_lev ) noexcept;
    amrex::MultiFab      * X ( int const a_lev ) noexcept;
    amrex::MultiFab      * tracer ( int const a_lev ) noexcept;

    amrex::MultiFab      * pert_p ( int const a_lev ) noexcept;
    amrex::MultiFab      * grad_p ( int const a_lev ) noexcept;

    amrex::MultiFab      * txfr ( int const a_lev ) noexcept;
    amrex::MultiFab      * chem_txfr ( int const a_lev ) noexcept;

    amrex::MultiFab      * vorticity ( int const a_lev ) noexcept;
    amrex::MultiFab      * div_tau ( int const a_lev ) noexcept;
    amrex::MultiFab      * mac_phi ( int const a_lev ) noexcept;

    amrex::MultiFab const* epf_old_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* rho_old_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* vel_old_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* T_old_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* h_old_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* X_old_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* tracer_old_const ( int const a_lev ) const noexcept;

    amrex::MultiFab const* epf_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* rho_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* vel_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* T_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* h_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* X_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* tracer_const ( int const a_lev ) const noexcept;

    amrex::MultiFab const* pert_p_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* grad_p_const ( int const a_lev ) const noexcept;

    amrex::MultiFab const* txfr_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* chem_txfr_const ( int const a_lev ) const noexcept;

    amrex::MultiFab const* vorticity_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* div_tau_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* mac_phi_const ( int const a_lev ) const noexcept;

    amrex::Array4<amrex::Real      > epf_old ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > rho_old ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > vel_old ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > T_old ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > h_old ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > X_old ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > tracer_old ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > epf ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > rho ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > vel ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > T ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > h ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > X ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > tracer ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > grad_p ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > pert_p ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > txfr ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > chem_txfr ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > vorticity ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > div_tau ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > mac_phi ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> epf_old_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> rho_old_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> vel_old_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> T_old_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> h_old_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> X_old_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> tracer_old_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> epf_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> rho_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> vel_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> T_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> h_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> X_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> tracer_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> pert_p_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> grad_p_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> txfr_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> chem_txfr_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> vorticity_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> div_tau_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> mac_phi_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Vector<std::unique_ptr<LevelData>> m_level_data;

  private:

    int m_finest_level = -1;

};

#endif
