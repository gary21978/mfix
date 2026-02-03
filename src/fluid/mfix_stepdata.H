
#ifndef _MFIX_STEPDATA_H_
#define _MFIX_STEPDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>

#include <mfix_bc.H>

class UniqueStepData {

  public:

    // Constructor
    UniqueStepData ( amrex::BoxArray                     const& a_ba,
                     amrex::DistributionMapping          const& a_dm,
                     amrex::FabFactory<amrex::FArrayBox> const& a_eb,
                     int const a_enthalpy,
                     int const a_nspecies,
                     int const a_ntracers,
                     int const a_ncomp_adv );

    // Destructor
    ~UniqueStepData () = default;

    std::unique_ptr<amrex::MultiFab> m_conv_velocity;
    std::unique_ptr<amrex::MultiFab> m_conv_scalars;
    std::unique_ptr<amrex::MultiFab> m_conv_species;

    std::unique_ptr<amrex::MultiFab> m_lap_T;
    std::unique_ptr<amrex::MultiFab> m_lap_tracer;

    std::unique_ptr<amrex::MultiFab> m_div_J;
    std::unique_ptr<amrex::MultiFab> m_div_hJ;

};

class SharedStepData {

  public:

    // Constructor
    SharedStepData ( amrex::BoxArray                     const& a_ba,
                     amrex::DistributionMapping          const& a_dm,
                     amrex::FabFactory<amrex::FArrayBox> const& a_eb,
                     int const a_enthalpy,
                     int const a_nspecies,
                     int const a_ntracers,
                     int const a_ncomp_adv,
                     int const a_nghost_state,
                     int const a_nghost_force,
                     int const a_nghost_mac );

    // Destructor
    ~SharedStepData () = default;

    std::unique_ptr<amrex::MultiFab> m_ep_u_mac;
    std::unique_ptr<amrex::MultiFab> m_ep_v_mac;
    std::unique_ptr<amrex::MultiFab> m_ep_w_mac;

    std::unique_ptr<amrex::MultiFab> m_RHS_proj;

    std::unique_ptr<amrex::MultiFab> m_depdt;

    std::unique_ptr<amrex::MultiFab> m_epf_nph;
    std::unique_ptr<amrex::MultiFab> m_rho_nph;

    std::unique_ptr<amrex::MultiFab> m_vel_S_p;
    std::unique_ptr<amrex::MultiFab> m_vel_S_c;

    std::unique_ptr<amrex::MultiFab> m_vel_forces;
    std::unique_ptr<amrex::MultiFab> m_tracer_forces;

    amrex::Array<std::unique_ptr<amrex::MultiFab>,3> m_J_k;
    amrex::Array<std::unique_ptr<amrex::MultiFab>,3> m_h_k;

    std::unique_ptr<amrex::MultiFab> m_eb_flow_velocity;
    std::unique_ptr<amrex::MultiFab> m_eb_flow_scalars;
    std::unique_ptr<amrex::MultiFab> m_eb_flow_species;

    std::unique_ptr<amrex::MultiFab> m_eb_temperature;

    std::unique_ptr<amrex::MultiFab> m_pert_pressure;

    std::unique_ptr<amrex::MultiFab> m_acceleration;

    std::unique_ptr<amrex::MultiFab> m_avg_particle_data;
};


class StepData {

  public:

    // Constructor
    StepData () = default;

    // Destructor
    ~StepData () = default;

    void define ( int const a_nlev,
                  int const a_enthalpy,
                  int const a_nspecies,
                  int const a_ntracers);

    int has_enthalpy () const noexcept { return m_has_enthalpy; }
    int has_tracer () const noexcept { return m_has_tracer; }
    int has_species () const noexcept { return m_has_species; }

    amrex::Real RHS_therm_p = 0.;

    amrex::Vector<std::unique_ptr<UniqueStepData>> m_step_data;

    // Methods that return vectors of constant MultiFabs

    amrex::Vector< amrex::MultiFab const* > conv_velocity_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > conv_scalars_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > conv_species_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > lap_T_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > lap_tracer_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > div_J_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > div_hJ_const () const noexcept;

    // Methods that return vectors of MultiFabs

    amrex::Vector< amrex::MultiFab      * > conv_velocity () const noexcept;
    amrex::Vector< amrex::MultiFab      * > conv_scalars () const noexcept;
    amrex::Vector< amrex::MultiFab      * > conv_species () const noexcept;

    amrex::Vector< amrex::MultiFab      * > lap_T () const noexcept;
    amrex::Vector< amrex::MultiFab      * > lap_tracer () const noexcept;

    amrex::Vector< amrex::MultiFab      * > div_J () const noexcept;
    amrex::Vector< amrex::MultiFab      * > div_hJ () const noexcept;

    // Methods that return constant Array4 pointers

    amrex::Array4<amrex::Real const> conv_velocity_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> conv_scalars_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> conv_species_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> lap_T_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> lap_tracer_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> div_J_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> div_hJ_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    // Methods that return Array4 pointers

    amrex::Array4<amrex::Real      > conv_velocity ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > conv_scalars ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > conv_species ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > lap_T ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > lap_tracer ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > div_J ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > div_hJ ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

  private:

    int m_nlev = 0;

    int m_has_enthalpy = 0;
    int m_has_species  = 0;
    int m_has_tracer   = 0;
};

class MFIXStepData {

  public:

    // Constructor
    MFIXStepData ( int const a_nlev,
                   amrex::Vector<amrex::BoxArray>                  const& a_ba,
                   amrex::Vector<amrex::DistributionMapping>       const& a_dm,
                   amrex::Vector<amrex::EBFArrayBoxFactory const*> const& a_factory,
                   int const a_nghost_state, int const a_nghost_force,
                   MFIXBoundaryConditions& a_bcs);

    // Destructor
    ~MFIXStepData () = default;


    void define ( amrex::Real const a_time,
                  int const a_density,  int const a_enthalpy,
                  int const a_nspecies, int const a_ntracers,
                  int const a_eb_flow,  int const a_eb_temperature);

    void include_corrector ();
    void include_pressure ();
    void include_acceleration ();
    void include_avg_particle_data (int ncomp);

    void set_pressure ( amrex::Vector< amrex::MultiFab const* > a_pressure );

    int has_eb_flow () const { return m_has_eb_flow; }
    int has_corrector () const { return m_has_corrector; }
    int has_eb_temperature () const { return m_has_eb_temperature; }
    int has_pressure () const { return m_has_pert_pressure; }
    int has_acceleration () const { return m_has_acceleration; }
    int has_avg_particle_data () const { return m_has_avg_particle_data; }

    StepData const& predictor_const () const noexcept { return m_predictor; }
    StepData const& corrector_const () const noexcept { return m_corrector; }

    StepData& predictor () noexcept { return m_predictor; }
    StepData& corrector () noexcept { return m_corrector; }

    // Methods that return vectors of constant MultiFabs

    amrex::Vector< amrex::Array<amrex::MultiFab const*,3>> mac_vel_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > ep_u_mac_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > ep_v_mac_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > ep_w_mac_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > RHS_proj_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > depdt_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > epf_nph_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > rho_nph_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > vel_S_p_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > vel_S_c_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > vel_forces_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > tracer_forces_const () const noexcept;

    amrex::Vector< amrex::Array<amrex::MultiFab const*,3> > J_k_const () const noexcept;
    amrex::Vector< amrex::Array<amrex::MultiFab const*,3> > h_k_const () const noexcept;

    amrex::Vector< amrex::MultiFab const*> eb_flow_velocity_const () const noexcept;
    amrex::Vector< amrex::MultiFab const*> eb_flow_scalars_const () const noexcept;
    amrex::Vector< amrex::MultiFab const*> eb_flow_species_const () const noexcept;
    amrex::Vector< amrex::MultiFab const*> eb_temperature_const () const noexcept;

    amrex::Vector< amrex::MultiFab const*> pert_pressure_const () const noexcept;
    amrex::Vector< amrex::MultiFab const*> acceleration_const () const noexcept;
    amrex::Vector< amrex::MultiFab const*> avg_particle_data_const () const noexcept;

    // Methods that return vectors of MultiFabs

    amrex::Vector< amrex::Array<amrex::MultiFab      *,3>> mac_vel () const noexcept;

    amrex::Vector< amrex::MultiFab      * > ep_u_mac () const noexcept;
    amrex::Vector< amrex::MultiFab      * > ep_v_mac () const noexcept;
    amrex::Vector< amrex::MultiFab      * > ep_w_mac () const noexcept;

    amrex::Vector< amrex::MultiFab      * > RHS_proj () const noexcept;

    amrex::Vector< amrex::MultiFab      * > depdt () const noexcept;

    amrex::Vector< amrex::MultiFab      * > epf_nph () const noexcept;
    amrex::Vector< amrex::MultiFab      * > rho_nph () const noexcept;

    amrex::Vector< amrex::MultiFab      * > vel_S_p () const noexcept;
    amrex::Vector< amrex::MultiFab      * > vel_S_c () const noexcept;

    amrex::Vector< amrex::MultiFab      * > vel_forces () const noexcept;
    amrex::Vector< amrex::MultiFab      * > tracer_forces () const noexcept;

    amrex::Vector< amrex::Array<amrex::MultiFab      *,3> > J_k () const noexcept;
    amrex::Vector< amrex::Array<amrex::MultiFab      *,3> > h_k () const noexcept;

    amrex::Vector< amrex::MultiFab      *> eb_flow_velocity () const noexcept;
    amrex::Vector< amrex::MultiFab      *> eb_flow_scalars () const noexcept;
    amrex::Vector< amrex::MultiFab      *> eb_flow_species () const noexcept;
    amrex::Vector< amrex::MultiFab      *> eb_temperature () const noexcept;

    amrex::Vector< amrex::MultiFab      *> pert_pressure () const noexcept;
    amrex::Vector< amrex::MultiFab      *> acceleration () const noexcept;
    amrex::Vector< amrex::MultiFab      *> avg_particle_data () const noexcept;

    // Methods that return constant MultiFab pointers

    amrex::Array<amrex::MultiFab const*,3> mac_vel_const ( int const a_lev ) const noexcept;

    amrex::MultiFab const* depdt_const ( int const a_lev ) const noexcept;

    amrex::MultiFab const* vel_S_p_const ( int const a_lev ) const noexcept;
    amrex::MultiFab const* vel_S_c_const ( int const a_lev ) const noexcept;

    // Methods that return MultiFab pointers

    amrex::Array<amrex::MultiFab      *,3> mac_vel ( int const a_lev ) const noexcept;

    amrex::MultiFab      * depdt ( int const a_lev ) const noexcept;

    amrex::MultiFab      * epf_nph ( int const a_lev ) const noexcept;
    amrex::MultiFab      * rho_nph ( int const a_lev ) const noexcept;

    // Methods that return constant Array4 pointers

    amrex::Array4<amrex::Real const> depdt_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> epf_nph_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> rho_nph_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> vel_S_p_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> vel_S_c_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> vel_forces_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> tracer_forces_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real const> acceleration_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real const> avg_particle_data_const ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    // Methods that return Array4 pointers

    amrex::Array4<amrex::Real      > depdt ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > epf_nph ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > rho_nph ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > vel_S_p ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > vel_S_c ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > vel_forces ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > tracer_forces ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    amrex::Array4<amrex::Real      > acceleration ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;
    amrex::Array4<amrex::Real      > avg_particle_data ( int const a_lev,
        const amrex::MFIter& a_mfi, int const a_scomp = 0 ) const noexcept;

    StepData m_predictor;

    StepData m_corrector;

    amrex::Vector<std::unique_ptr<SharedStepData>> m_shared_data;

  private:

    int const m_nlev;

    amrex::Vector<amrex::BoxArray>            const& m_ba;
    amrex::Vector<amrex::DistributionMapping> const& m_dm;
    amrex::Vector<amrex::EBFArrayBoxFactory const*> m_factory;

    int const m_nghost_state;
    int const m_nghost_force;
    int const m_nghost_mac = 1;

    // Boundary conditions
    MFIXBoundaryConditions& m_bcs;

    int m_defined = 0;

    int m_density  = 0;
    int m_enthalpy = 0;
    int m_nspecies = 0;
    int m_ntracers = 0;

    int m_ncomp_adv = 0;

    int m_has_eb_flow        = 0;
    int m_has_corrector      = 0;
    int m_has_eb_temperature = 0;
    int m_has_pert_pressure  = 0;

    int m_has_acceleration = 0;

    int m_has_avg_particle_data = 0;

    void include_eb_flow ( amrex::Real const a_time );
    void include_eb_temperature ( amrex::Real const a_time );

};

#endif
