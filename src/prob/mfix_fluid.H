#ifndef _MFIX_FLUID_H_
#define _MFIX_FLUID_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>

#include <mfix_fluid_var.H>
#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_properties.H>


// Forward declaration
class MFIXFluidPhase;

// Advection Type
enum struct AdvectionType { Undefined=-1, MOL, Godunov };

struct AdvectionComp { enum { density = 0, enthalpy, tracer }; };

struct InterphaseTxfrIndexes {

  InterphaseTxfrIndexes ()
    : drag_coeff(0)
    , vm_coeff(1)
    , vel_src_c(2)
    , convection_coeff(5)
    , energy_source(6)
    , count(7)
  {}

  int const drag_coeff;
  int const vm_coeff;
  int const vel_src_c;
  int const convection_coeff;
  int const energy_source;
  int const count;
};


/**
 * @struct DualGridAuxIndexes
 * @brief The start index for each component used to transfer grid data
 * to particle positions.
 *
 * A local MultiFab is created when grid values are mapped to particle
 * particle positions with grid data copied to specific components
 * of the destination MultiFab. For example, fluid velocity is copied
 * to components (0,1,2) and fluid volume fraction is copied to
 * component 4. Components are omitted or included based on simulation
 * settings.
 */
struct DualGridAuxIndexes {
  DualGridAuxIndexes ( const int a_include_energy
                     , const int a_nspecies
                     , const int a_is_mixture )
    : include_energy(a_include_energy)
    , nspecies(a_nspecies)
    , is_mixture(a_is_mixture)
    , vel_g(0)
    , ep_g(vel_g+3)
    , ro_g(ep_g+1)
    , T_g(ro_g+1)
    , h_g(T_g+1*include_energy)
    , X_gk(h_g+1*include_energy)
    , count(3 + 1 + 1 + (1+1)*include_energy + nspecies)
  {}

  int const include_energy;
  int const nspecies;
  int const is_mixture;

  int const vel_g;
  int const ep_g;
  int const ro_g;
  int const T_g;
  int const h_g;
  int const X_gk;
  int const count;
};

enum struct EddyViscModel {None=0, Smagorinsky, WALE, Usr};
enum struct SuspViscModel {None=0, Einstein, Brinkman, Roscoe, ChengLaw, Sato, Subramaniam};

struct MFIXFluidConstraint {

    int isIncompressibleFluid () const
    { return (m_constraint == m_IncompressibleFluid); }

    int isIdealGasOpenSystem () const
    { return (m_constraint == m_IdealGasOpenSystem); }

    int isIdealGasClosedSystem () const
    { return (m_constraint == m_IdealGasClosedSystem); }

    int isIdealGas () const
    { return (isIdealGasOpenSystem() || isIdealGasClosedSystem()); }

    int include_depdt () const { return m_include_depdt; }

    int set ( std::string const& a_constraint ) {

      std::string constraint = amrex::toLower(a_constraint);

      if (constraint == "incompressiblefluid") {
        m_constraint = m_IncompressibleFluid;

      } else if (constraint == "idealgasopensystem") {
        m_constraint = m_IdealGasOpenSystem;

      } else if (constraint == "idealgasclosedsystem") {
        m_constraint = m_IdealGasClosedSystem;

      } else { // Invalid. Flag error.
        return 1;
      }
      return 0;
    }

    void set_depdt ( int const a_flag ) { m_include_depdt = a_flag; }

  private:

    enum { m_Undefined,
           m_IncompressibleFluid,
           m_IdealGasOpenSystem,
           m_IdealGasClosedSystem
         };

    int m_constraint = m_Undefined;

    int m_include_depdt = 0;

};


class MFIXFluidParms {

  public:
    //

    AMREX_GPU_HOST_DEVICE
    MFIXFluidParms ();

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int nspecies () const { return m_nspecies; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_MW_gk (const int n_g) const
    {
      AMREX_ASSERT((n_g >= 0 && n_g < m_nspecies) || (n_g == 0 && m_nspecies == 0));
      return m_MW_gk[n_g];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_MW_g () const { return get_MW_gk(0); }

    // Compute fluid thermal conductivity
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_k_g (const amrex::Real /*T_g*/) const { return m_k_g; }

    // Compute fluid thermal conductivity
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_partial_k_g (const amrex::Real /*T_g*/) const { return 0; }

    // Compute fluid species diffusivity
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_D_g () const { return m_D_g; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_lagrangian_stoich_coeff (const int n_g,
                                          const int q) const
    {
      AMREX_ASSERT(m_nspecies > 0);
      AMREX_ASSERT(n_g >= 0 && n_g < m_nspecies);
      AMREX_ASSERT(m_lagrangian_nreactions > 0);
      AMREX_ASSERT(q >= 0 && q < m_lagrangian_nreactions);
      return m_lagrangian_stoich_coeffs[n_g*m_lagrangian_nreactions+q];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_eulerian_stoich_coeff (const int n_g,
                                        const int q) const
    {
      AMREX_ASSERT(m_nspecies > 0);
      AMREX_ASSERT(n_g >= 0 && n_g < m_nspecies);
      AMREX_ASSERT(m_eulerian_nreactions > 0);
      AMREX_ASSERT(q >= 0 && q < m_eulerian_nreactions);
      return m_eulerian_stoich_coeffs[n_g*m_eulerian_nreactions+q];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int species_idx (const char* species_name) const
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXFluidParms::species_idx()");
#endif

      int ordered_idx = species_algorithms::linear_search(m_species_names,
          m_species_names_idxs, m_nspecies, species_name);
      return m_species_names_IDs[ordered_idx];
    }

  private:
    // Friendship so MFIXFluidPhase can access private members
    friend class MFIXFluidPhase;

    amrex::Real  m_k_g;
    int          m_nspecies;
    char*        m_species_names;
    int*         m_species_names_idxs;
    int*         m_species_names_IDs;
    amrex::Real* m_MW_gk;
    amrex::Real  m_D_g;
    int          m_lagrangian_nreactions;
    int          m_eulerian_nreactions;
    amrex::Real* m_lagrangian_stoich_coeffs;
    amrex::Real* m_eulerian_stoich_coeffs;
};


class MFIXFluidPhase {

  public:
    //
    struct ThermalConductivityModel { enum {Undefined=-1, Constant}; };


    // Ideal gas constant
    static constexpr amrex::Real R = 8.31446261815324;

    // Constructor
    MFIXFluidPhase();

    // Destructor
    ~MFIXFluidPhase();

    int Initialize (const MFIXSpecies& species, const MFIXReactions& reactions);

    const std::string& names (const int i) const { return m_names[i]; }

    bool solve () const { return (m_solve == 1);  }
    bool solve_density () const { return (m_solve_density == 1);  }
    bool solve_tracer () const { return (m_solve_tracer == 1);  }
    bool solve_enthalpy () const { return (m_solve_enthalpy == 1);  }
    bool solve_species () const { return (m_solve_species == 1);  }

    const amrex::Vector<std::string>& species_names () const { return m_species_names; }
    const std::string& species_names (const int n_g) const { return m_species_names[n_g]; }
    int nspecies () const { return m_nspecies; }

    int ntracer () const { return m_ntracer; }

    amrex::Vector<amrex::Real> tracer_diff_coeff () const
    { return m_tracer_diff_coeff; }

    int isInitialized () const { return m_is_initialized; }
    int isMixture () const { return m_is_a_mixture; }

    template <amrex::RunOn run_on>
    const MFIXFluidParms& parameters () const
    {
      if constexpr (run_on == amrex::RunOn::Cpu)
        return *m_h_parameters;
      else if constexpr (run_on == amrex::RunOn::Gpu)
        return *m_d_parameters;
    }

    ThermoProperties props;

    MFIXFluidConstraint constraint;

    amrex::Real tracer () const { return m_trac_0; }
    int p_therm_defined () const { return m_p_therm_defined; }
    amrex::Real thermodynamic_pressure () const { return m_thermodynamic_pressure; }

    void set_p_therm_defined (const int i) { m_p_therm_defined = i; }
    void set_thermodynamic_pressure (const amrex::Real p) { m_thermodynamic_pressure = p; }

    EddyViscModel eddy_visc_model () const { return m_eddy_visc_model; }

    bool eddy_visc_model ( const EddyViscModel a_model ) const
    { return (a_model == m_eddy_visc_model); }

    SuspViscModel susp_visc_model () const { return m_susp_visc_model; }

    bool susp_visc_model ( SuspViscModel a_model ) const
    { return (a_model == m_susp_visc_model); }

    amrex::Real smagorinsky_constant () const { return m_smagorinsky_constant; }
    amrex::Real wale_constant () const { return m_wale_constant; }
    amrex::Real brinkman_constant () const { return m_brinkman_constant; }
    amrex::Real roscoe_c1 () const { return m_roscoe_c1; }
    amrex::Real roscoe_c2 () const { return m_roscoe_c2; }
    amrex::Real chenglaw_constant () const { return m_chenglaw_constant; }
    amrex::Real sato_constant () const { return m_sato_constant; }
    amrex::Real subramaniam_constant () const { return m_subramaniam_constant; }
    amrex::Real max_effective_viscosity_factor () const { return m_max_effective_viscosity_factor; }

  private:

    int m_ntypes;

    // List of recognized fluid viscosity models
    EddyViscModel m_eddy_visc_model = EddyViscModel::None;
    SuspViscModel m_susp_visc_model = SuspViscModel::None;

    // List of recognized fluid thermal conductivity models
    int m_thermal_conductivity_model = ThermalConductivityModel::Undefined;

    // Names of fluids
    amrex::Vector<std::string> m_names;

    // Flag to solve fluid equations
    int m_solve;

    // Flag to solve fluid density equation
    int m_solve_density;

    // Flag to solve fluid tracer equation
    int m_solve_tracer;

    // Number of tracers
    int m_ntracer;

    // Tracer diffusion coefficients
    amrex::Vector<amrex::Real> m_tracer_diff_coeff;

    // Specified constant tracer value
    amrex::Real m_trac_0;

    // Flag to solve fluid enthalpy equation
    int m_solve_enthalpy;

    // Specified constant gas phase thermal conductivity coefficient
    amrex::Real m_k_g0;

    // Specified gas phase thermodynamic pressure
    amrex::Real m_thermodynamic_pressure;
    int m_p_therm_defined;

    // Flag to solve species fluid equations
    int m_solve_species;

    // Smagorinsky-Lilly constant
    amrex::Real m_smagorinsky_constant;

    // WALE constant
    amrex::Real m_wale_constant;

    // Brinkman constant
    amrex::Real m_brinkman_constant;

    // Roscoe constants
    amrex::Real m_roscoe_c1;
    amrex::Real m_roscoe_c2;

    // ChengLaw constant
    amrex::Real m_chenglaw_constant;

    // Sato constant
    amrex::Real m_sato_constant;

    // Subramaniam constant
    amrex::Real m_subramaniam_constant;

    // Factor to limit max of the effective viscosity
    amrex::Real m_max_effective_viscosity_factor;

    // Fluid phase species names
    amrex::Vector<std::string> m_species_names;
    amrex::Gpu::HostVector<char> m_h_species_names;
    amrex::Gpu::HostVector<int> m_h_species_names_idxs;
    amrex::Gpu::HostVector<int> m_h_species_names_IDs;
    amrex::Gpu::DeviceVector<char> m_d_species_names;
    amrex::Gpu::DeviceVector<int> m_d_species_names_idxs;
    amrex::Gpu::DeviceVector<int> m_d_species_names_IDs;

    // Species unique identifying code
    amrex::Vector<int> m_species_IDs;

    // Total number of fluid species
    int m_nspecies;

    // Specified constant gas phase species molecular weight
    amrex::Gpu::HostVector<amrex::Real> m_MW_gk0;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_MW_gk0;

    // Specified constant gas phase species diffusion coefficients
    amrex::Real m_D_g0;

    // Flag to understand if fluid is a mixture of fluid species
    int m_is_a_mixture;

    // Tables for stoichiometric coefficients for each species in each
    // lagrangian reaction
    amrex::Gpu::HostVector<amrex::Real> m_lagrangian_stoich_coeffs;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_lagrangian_stoich_coeffs;

    // Tables for stoichiometric coefficients for each species in each
    // eulerian reaction
    amrex::Gpu::HostVector<amrex::Real> m_eulerian_stoich_coeffs;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_eulerian_stoich_coeffs;

    // Fluid parameters for GPU access
    MFIXFluidParms* m_h_parameters;
    MFIXFluidParms* m_d_parameters;

    int m_is_initialized;

}; // class MFIXFluidPhase


struct FLUID_t {

  FLUID_t(MFIXFluidPhase* fluid_in)
    : fluid(fluid_in)
    , flow_thru_eb(false)
    , volfrac(-1.0)
    , pressure(-1.0)
    , pressure_defined(0)
    , density({1})
    , temperature({1})
    , massflow({1})
    , volflow({1})
    , velocity({1,3})
  { }

  MFIXFluidPhase* fluid;

  bool flow_thru_eb;

  amrex::Real volfrac;
  amrex::Real pressure;
  int pressure_defined; // Flag to check if pressure is defined

  FVAR_ density;
  FVAR_ temperature;
  amrex::Vector<FVAR_> species;
  FVAR_ massflow;
  FVAR_ volflow;
  FVAR_ velocity;

  amrex::Real get_density (amrex::Real a_time = 0) const;

  amrex::Real get_temperature (amrex::Real a_time = 0) const;

  amrex::Real get_species (int a_n, amrex::Real a_time = 0) const;

  amrex::Real get_massflow (amrex::Real a_time = 0) const;

  amrex::Real get_volflow (amrex::Real a_time = 0) const;

  amrex::Real get_velocity_mag (amrex::Real time = 0) const;

  amrex::RealVect get_velocity (amrex::Real a_time = 0) const;

};


#endif
