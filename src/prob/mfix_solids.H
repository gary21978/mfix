#ifndef MFIX_SOLIDS_PARMS_H_
#define MFIX_SOLIDS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_RandomEngine.H>

#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_distributions.H>
#include <mfix_properties.H>

// Forward declaration
class MFIXSolidsParms;


struct SoArealData
{
  enum {
    radius = 0,
    density,
    velx,
    vely,
    velz,
    omegax,
    omegay,
    omegaz,
    statwt,
    drag_coeff,
    vel_source_x,
    vel_source_y,
    vel_source_z,
    temperature,
    count
  };

  static
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real volume (const amrex::Real radius)
  {
    constexpr amrex::Real PI_by_4_over_3 = (4./3.) * M_PI;
    return PI_by_4_over_3 * (radius*radius*radius);
  }

  static
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real oneOverI (const amrex::Real radius,
                        const amrex::Real mass)
  { return 2.5 / (mass * (radius*radius)); }
};

/**
 * @brief Class to identify real data stored at runtime
 *
 * This structure contains indexes for particle quantities that are
 * model dependent. The data is stored as a struct of arrays (SoA).
 */
class runtimeRealData
{
  public:

    /**
     * @brief Default constructor for when there are no solids.
     */
    runtimeRealData () = default;

    /**
     * @brief Constructor to initialize particle runtime reals layout.
     *
     * @param a_nspecies Number of species comprising particles
     * @param a_pic_solve Flag for PIC model; include volume fraction, ep_s
     * @param a_include_vm Flag to include virtual mass coefficient
     * @param a_include_acc Flag to include particle acceleration
     * @param a_include_energy Flag to include energy variables
     * @param a_max_contacts Max number of simultaneous contacts
     */
    AMREX_GPU_HOST_DEVICE
    runtimeRealData (int const a_nspecies,
                     int const a_pic_solve,
                     int const a_include_vm,
                     int const a_include_acc,
                     int const a_include_energy,
                     int const a_max_contacts)
      : ncomp_species((a_nspecies>0) ? a_nspecies : 0)
      , ncomp_eps(a_pic_solve         ? 1 : 0)
      , ncomp_vm(a_include_vm         ? 1 : 0)
      , ncomp_acc(a_include_acc       ? 3 : 0)
      , ncomp_energy(a_include_energy ? 2 : 0)
      , max_contacts(a_max_contacts)
      // Set indices of runtime real objects
      , X_sn(0)
      , ep_s(X_sn + ncomp_species)
      , vm_coeff(ep_s + ncomp_eps)
      , acceleration(vm_coeff + ncomp_vm)
      , energy_source(acceleration + ncomp_acc)
      , conv_coeff(acceleration + ncomp_acc + 1)
      , pft_neighbor_idx(energy_source + ncomp_energy)
      // Number of entries in the runtime real data
      , count_no_tan_history(energy_source + ncomp_energy)
      , count(count_no_tan_history + 3*a_max_contacts)
    {}

    const int ncomp_species = 0; ///< Number of species comprising particles
    const int ncomp_eps     = 0; ///< Number of volume fraction components
    const int ncomp_vm      = 0; ///< Number of virtual mass components
    const int ncomp_acc     = 0; ///< Number of acceleration components
    const int ncomp_energy  = 0; ///< Number of energy components

    const int max_contacts  = 0; ///< Max number of simultaneous contacts

    const int X_sn             = -1; ///< Index of first species (start component)
    const int ep_s             = -1; ///< Index of volume fraction
    const int vm_coeff         = -1; ///< Index of virtual mass coefficient
    const int acceleration     = -1; ///< Index of first acceleration component
    const int energy_source    = -1; ///< Index of total, explicit energy source
    const int conv_coeff       = -1; ///< Index of convection coefficient
    const int pft_neighbor_idx = -1; ///< Index of first contact history component

    const int count_no_tan_history = 0; ///< total count excluding contact history
    const int count                = 0; ///< total count including contact history

    /// Identify if runtime real data contains species
    int contains_species () const noexcept { return (ncomp_species > 0) ? 1 : 0; }

    /// Identify if runtime real data contains volume fraction
    int contains_eps () const noexcept { return (ncomp_eps > 0) ? 1 : 0; }

    /// Identify if runtime real data contains virtual mass coefficient
    int contains_vm () const noexcept { return (ncomp_acc > 0) ? 1 : 0; }

    /// Identify if runtime real data contains particle acceleration
    int contains_acc () const noexcept { return (ncomp_vm > 0) ? 1 : 0; }

    /// Identify if runtime real data contains energy variables
    int contains_energy () const noexcept { return (ncomp_energy > 0) ? 1 : 0; }

};

//Particle int data stored at runtime in a SoA style
class runtimeIntData
{
   public:
    AMREX_GPU_HOST_DEVICE
    runtimeIntData (const int contacts_max_tan_history)
      : max_contacts_tan_history(contacts_max_tan_history)
      , cpu_id_idx(0)
      , pid_idx(cpu_id_idx + contacts_max_tan_history)
      , wall_pid_idx(-2) // special index for wall TODO: handle better
      , count(2*contacts_max_tan_history)
   {}

    const int max_contacts_tan_history;
    const int cpu_id_idx;
    const int pid_idx;
    const int wall_pid_idx;
    const int count;
};

struct SoAintData
{
  enum {
    phase = 0,
    state,
#if MFIX_POLYDISPERSE
    ptype,
#endif
    count
  };
};

struct state
{
  enum { invalid=-1, entering=0, normal=1 };
};


/**
 * \class MFIXSolidsPhase
 * \brief Settings defining the solids
 */
class MFIXSolidsPhase {

  public:

    /// Simple constructor that initializes member variables
    MFIXSolidsPhase ();

    /// Destructor used to deallocate MFIXSolidsParms pointers
    ~MFIXSolidsPhase ();

    /**
     * Available models for computing solids thermal conductivity.
     * \param Undefined used to initialize to invalid state.
     * \param Constant constant thermal conductivity.
     */
    struct ThermalConductivityModel { enum {Undefined=-1, Constant}; };
    struct DensityModel { enum {Undefined=-1, SpeciesMixture}; };

    // TODO: remove phase <--> index conversion as this is a very old
    // artifact when moving from FORTRAN to C++
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static int phase_to_index (const int phase)
    { return phase-1; }

    void Initialize (const MFIXSpecies& species,
                     const MFIXReactions& reactions);

    /// Test if the class was previously initialized
    int isInitialized () const { return m_is_initialized; }

    int ntypes () const { return m_ntypes; }

    const amrex::Vector<std::string>& names () const { return m_names; }
    const std::string& names (const int n) const { return m_names[n]; }

    const amrex::Vector<int>& phases () const { return m_phases; }
    int phases (const int n) const { return m_phases[n]; }

    /// Return phase ID given a solids name. This is used in processing inputs data
    int name_to_phase (const std::string& name) const;

    /// Flag indicating species equations are being solved.
    int solve_species () const { return m_solve_species; }

    /// Flag indicating energy equations are being solved.
    int solve_enthalpy () const { return m_solve_enthalpy; }

    /// User defined, constant energy source for solids
    amrex::Real enthalpy_source () const { return m_enthalpy_source; }

    /// Flag indicating that particle positions and velocity are to be updated.
    int update_momentum () const { return m_update_momentum; }

    /// Flag indicating that particle mass is to be updated.
    int update_mass () const { return m_update_mass; }

    /// Flag indicating that particle temperature is to be updated.
    int update_enthalpy () const { return m_update_enthalpy; }

    /// number of species comprising the solids
    int nspecies () const { return m_nspecies; }

    /// returns true if the solids is comprised of more than one species
    int isMixture () const { return m_is_a_mixture; }

    /// returns true if particle acceleration is stored as runtime variable
    int track_acceleration () const { return m_track_acceleration; }

    const amrex::Vector<std::string>& species_names () const { return m_species_names; }
    const std::string& species_names (const int n) const { return m_species_names[n]; }

    const amrex::Gpu::HostVector<amrex::Real>& species_densities () const { return m_h_species_densities; }
    amrex::Real species_density (const int n) const { return m_h_species_densities[n]; }
    short density_model () const { return m_density_model; }

    template <amrex::RunOn run_on>
    const MFIXSolidsParms& parameters () const
    { return (run_on == amrex::RunOn::Cpu) ? *m_h_parameters : *m_d_parameters; }

    ThermoProperties props;

  private:

    int m_ntypes;

    /// List of recognized solids thermal conductivity models
    int m_thermal_conductivity_model;

    /// List of recognized solids species densities models
    int m_density_model = DensityModel::Undefined;

    /// Names of solids used in IC/BC setups
    amrex::Vector<std::string> m_names;

    /// Phases IDs for each solids type
    amrex::Vector<int> m_phases;

    /// Flag to turn on/off particles mass update
    int m_update_mass;

    /// Flag to solve species
    int m_solve_species;

    /// Flag to turn on/off particles velocity update
    int m_update_momentum;

    /// User-defined enthalpy volumetric source for particles
    amrex::Real m_enthalpy_source;

    /// Flag to turn on/off particles enthalpy update
    int m_update_enthalpy;

    /// Flag to solve enthalpy
    int m_solve_enthalpy;

    /// Particle species
    amrex::Vector<std::string> m_species_names;
    amrex::Gpu::HostVector<char> m_h_species_names;
    amrex::Gpu::HostVector<int> m_h_species_names_idxs;
    amrex::Gpu::HostVector<int> m_h_species_names_IDs;
    amrex::Gpu::DeviceVector<char> m_d_species_names;
    amrex::Gpu::DeviceVector<int> m_d_species_names_idxs;
    amrex::Gpu::DeviceVector<int> m_d_species_names_IDs;

    /// Species unique identifying code
    amrex::Vector<int> m_species_IDs;

    /// Total number of dem species
    int m_nspecies;

    // Specified constant solid phase species molecular weight
    amrex::Gpu::HostVector<amrex::Real> m_MW_sn0;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_MW_sn0;

    // Specified constant solid phase species densities
    amrex::Gpu::HostVector<amrex::Real> m_h_species_densities = {};
    amrex::Gpu::DeviceVector<amrex::Real> m_d_species_densities = {};

    // Flag to set solids as a mixture of solids species
    int m_is_a_mixture;

    // Flag indicating that particle acceleration is stored in the particle container.
    int m_track_acceleration;

    // Tables for stoichiometric coefficients for each species in each reaction
    amrex::Gpu::HostVector<amrex::Real> m_lagrangian_stoich_coeffs;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_lagrangian_stoich_coeffs;

    // Specified conductivity
    amrex::Gpu::HostVector<amrex::Real> m_kp_sn0;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_kp_sn0;

    // Specified fluid-lens constant
    amrex::Real m_flpc;

    // Specified roughness constant
    amrex::Real m_rough;

    // Flag to do pfp conduction
    int m_do_conduction;

    // Solids parameters for GPU access
    MFIXSolidsParms* m_h_parameters;
    MFIXSolidsParms* m_d_parameters;

    int m_is_initialized;
};


class MFIXSolidsParms {

  public:
    AMREX_GPU_HOST_DEVICE
    MFIXSolidsParms ()
      : m_nspecies(0)
      , m_solve_dem(1)
      , m_species_names(nullptr)
      , m_species_names_idxs(nullptr)
      , m_species_names_IDs(nullptr)
      , m_MW_sn(nullptr)
      , m_lagrangian_nreactions(0)
      , m_lagrangian_stoich_coeffs(nullptr)
      , m_kp_sn(nullptr)
      , m_flpc(0.4)
      , m_rough(2E-8)
      , m_do_pfp_cond(0)
      , m_thermal_conductivity_model(MFIXSolidsPhase::ThermalConductivityModel::Undefined)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int nspecies () const { return m_nspecies; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int solve_dem () const { return m_solve_dem; }

    // Driver molecular weight
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_MW_s () const { return get_MW_sn(0); }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_MW_sn (const int n_s) const
    {
      AMREX_ASSERT((n_s >= 0 && n_s < m_nspecies) || (n_s == 0 && m_nspecies == 0));
      return m_MW_sn[n_s];
    }

    // Driver solid density
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_density (const int n_s) const
    {
      AMREX_ASSERT((n_s >= 0 && n_s < m_nspecies) || (n_s == 0 && m_nspecies == 0));
      return m_species_densities[n_s];
    }

    // Get lagrangian reaction coefficients
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_lagrangian_stoich_coeff (const int n_s, const int q) const
    {
      AMREX_ASSERT(m_nspecies > 0);
      AMREX_ASSERT(n_s >= 0 && n_s < m_nspecies);
      AMREX_ASSERT(m_lagrangian_nreactions > 0);
      AMREX_ASSERT(q >= 0 && q < m_lagrangian_nreactions);
      return m_lagrangian_stoich_coeffs[n_s*m_lagrangian_nreactions+q];
    }

    // Get fluid lens proportionality constant
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_flpc () const { return m_flpc; }

    // Get min conduction distance
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_min_cond () const { return m_rough; }

    // Do we enter PFP conduction?
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int get_do_pfp_cond () const { return m_do_pfp_cond; }

    //************************************************************

    // Calc conductivity
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_kp_sn (const amrex::Real, const int n_s) const
    {
      AMREX_ASSERT((n_s >= 0 && n_s < m_nspecies) || (n_s == 0 && m_nspecies == 0));
      amrex::Real result(0);

      if (m_thermal_conductivity_model == MFIXSolidsPhase::ThermalConductivityModel::Constant) {
        result = m_kp_sn[n_s];
      }

      return result;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int species_idx (const char* species_name) const
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXSolidsParms::species_idx()");
#endif

      int ordered_idx = species_algorithms::linear_search(m_species_names,
          m_species_names_idxs, m_nspecies, species_name);
      return m_species_names_IDs[ordered_idx];
    }

  private:
    friend class MFIXSolidsPhase;

    int   m_nspecies;
    int   m_solve_dem;
    char* m_species_names;
    int*  m_species_names_idxs;
    int*  m_species_names_IDs;
    amrex::Real* m_MW_sn;
    amrex::Real* m_species_densities = nullptr;
    int  m_lagrangian_nreactions;
    amrex::Real* m_lagrangian_stoich_coeffs;
    amrex::Real* m_kp_sn;
    amrex::Real m_flpc;
    amrex::Real m_rough;
    int m_do_pfp_cond;
    int m_thermal_conductivity_model;
};


struct SOLIDS_t {

  public:

    SOLIDS_t() :
      phase(0),
      volfrac(0.0),
      points_at_pack(1.0),
      temperature(0.0),
      velocity(0),
      velmag(0),
      volflow(0.),
      vol_remainder(0.),
      species(0),
      diameter(),
      density(),
      m_diameter_data()
    {}

  std::string name;

  int phase;

  amrex::Real volfrac;
  amrex::Real points_at_pack;
  amrex::Real temperature;
  amrex::Vector<amrex::Real> velocity;
  amrex::Real velmag;
  amrex::Real volflow;
  amrex::Real vol_remainder;

  amrex::Vector<SPECIES_t> species;

  INPUT_DIST_t diameter;

  INPUT_DIST_t density;

  long calc_particles_from_vol (amrex::Real const a_vol,
                                amrex::Real const a_multi_particle_volume,
                                int const include_volfrac = 1) const;

  const INPUT_DIST_t* get_diameter () const { return &diameter; }

  void copyAsync () { diameter.copyAsync(&m_diameter_data); }

  INPUT_DIST_DATA_t* get_diameter_data () { return &m_diameter_data; }

  private:

    INPUT_DIST_DATA_t m_diameter_data;

};

#endif
