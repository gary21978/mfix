#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>

#include <mfix_reporter.H>
#include <mfix_solids.H>
#include <mfix_species.H>
#include <mfix_fluid.H>
#include <mfix_reporter.H>


using namespace amrex;


MFIXSolidsPhase::MFIXSolidsPhase()
  : m_ntypes(0)
  , m_thermal_conductivity_model(ThermalConductivityModel::Undefined)
  , m_names(0)
  , m_update_mass(1)
  , m_solve_species(0)
  , m_update_momentum(1)
  , m_enthalpy_source(0)
  , m_update_enthalpy(1)
  , m_solve_enthalpy(0)
  , m_species_names(0)
  , m_h_species_names(0)
  , m_h_species_names_idxs(0)
  , m_d_species_names(0)
  , m_d_species_names_idxs(0)
  , m_species_IDs(0)
  , m_nspecies(0)
  , m_MW_sn0(0)
  , m_d_MW_sn0(0)
  , m_is_a_mixture(0)
  , m_track_acceleration(0)
  , m_lagrangian_stoich_coeffs(0)
  , m_d_lagrangian_stoich_coeffs(0)
  , m_kp_sn0(0)
  , m_d_kp_sn0(0)
  , m_flpc(0.4)
  , m_rough(2.E-8)
  , m_do_conduction(0)
  , m_h_parameters(nullptr)
  , m_d_parameters(nullptr)
  , m_is_initialized(0)
{}


MFIXSolidsPhase::~MFIXSolidsPhase()
{
  if (m_h_parameters != nullptr) { delete m_h_parameters; }
  if (m_d_parameters != nullptr) { delete m_d_parameters; }
}


int
MFIXSolidsPhase::name_to_phase (const std::string& name) const
{
  if (m_ntypes == 0) {
    Print() << "Can't look for solid type. No solids types found\n";
    amrex::Abort("Error. Fix inputs");
  }

  for (int n(0); n < m_ntypes; ++n) {
    if (m_names[n].compare(name) == 0) {
      return m_phases[n];
    }
  }

  amrex::Abort("Solid name provided does not exist");

  return -1;
}


void
MFIXSolidsPhase::Initialize (const MFIXSpecies& species,
                             const MFIXReactions& reactions)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.isInitialized(),
      "Species not initialized. Can't initialize solids phase before species initialization");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(reactions.isInitialized(),
      "MFIXReactions not initialized. Can't initialize solids phase before reactions initialization");

  // Flag for initialization
  m_is_initialized = 1;

  amrex::ParmParse pp_solids("solids");

  pp_solids.queryarr("types", m_names);
  m_ntypes = m_names.size();

  // Disable the solids solver if the solids is defined as "" or "None"
  if ((m_ntypes > 0) && (amrex::toLower(m_names[0]).compare("none") != 0)) {

    // Resize phases and assign integer IDs to each phase as {1,2,...,ntypes+1}
    m_phases.resize(m_ntypes);
    for (int n(0); n < m_ntypes; ++n)
      m_phases[n] = n+1;

    // Get mfix global inputs ------------------------------------------------//
    amrex::ParmParse ppMFIX("mfix");

    ppMFIX.query("advect_enthalpy", m_solve_enthalpy);

    ppMFIX.query("solve_species", m_solve_species);

    // Query updating flags for particles mass, momentum and enthalpy
    {
      amrex::ParmParse pp_particles("mfix.particles");

      pp_particles.query("update_mass", m_update_mass);
      m_update_mass = m_update_mass && m_solve_species;

      pp_particles.query("update_momentum", m_update_momentum);

      pp_particles.query("enthalpy_source", m_enthalpy_source);

      pp_particles.query("update_enthalpy", m_update_enthalpy);
      m_update_enthalpy = m_update_enthalpy && m_solve_enthalpy;
    }

    if (m_solve_species) {
      // Query solids species
      pp_solids.queryarr("species", m_species_names);

      m_nspecies = m_species_names.size();

      if (m_nspecies > 0) {

        std::map<std::string, int> ordered_species_names;
        int length(0);
        for (int n(0); n < m_nspecies; ++n) {
          const std::string name = m_species_names[n];
          ordered_species_names.insert(std::make_pair(name, n));
          length += name.size();
        }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nspecies == Long(ordered_species_names.size()),
            "There are solids species defined multiple times. Check inputs");

        m_h_species_names.clear();
        m_h_species_names.reserve(length+1);

        m_h_species_names_idxs.clear();
        m_h_species_names_idxs.reserve(m_nspecies+1);

        m_h_species_names_IDs.clear();
        m_h_species_names_IDs.reserve(m_nspecies);

        length = 0;
        for (auto elem: ordered_species_names) {
          const std::string& name = get<0>(elem);
          for (auto character: name) {
            m_h_species_names.push_back(character);
          }
          m_h_species_names_idxs.push_back(length);
          m_h_species_names_IDs.push_back(get<1>(elem));
          length += name.size();
        }
        m_h_species_names.push_back('\0');
        m_h_species_names_idxs.push_back(length);
      }
    }

    if (!m_solve_species || m_species_names.size() == 0 ||
        amrex::toLower(m_species_names[0]).compare("none") == 0)
    {
      m_species_names.clear();
      m_h_species_names.clear();
      m_h_species_names_idxs.clear();
      m_solve_species = 0;
      m_nspecies = 0;
    } else {
      m_solve_species = 1;
      m_nspecies = m_species_names.size();
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_species_names.size() > 0,
          "No input provided for solids.species");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nspecies <= species.nspecies(),
          "Solids species number is higher than total species number");
    }

    // Flag to determine if we want to solve the solid as a mixture
    m_is_a_mixture = static_cast<int>(m_nspecies > 1);

    // Get specific heat and enthalpy inputs ---------------------------//
    if (m_solve_enthalpy) {
      if (props.specificHeat.define("solids") != 0)
      {
        Abort("Failed to initialize specific heat!");
      }
    }

    // Get solids species parameters from species class
    if (m_solve_species) {

      std::string density_model("None");
      pp_solids.query("density", density_model);

      m_species_IDs.resize(m_nspecies);
      m_MW_sn0.resize(m_nspecies);

      if (amrex::toLower(density_model) == "mixture") {

        if (species.densities_k().size() == 0) {
          reporter::Log(reporter::Error, __FILE__, __LINE__)
            << "Solids density model specified as mixture\n"
            << "but species constant densities have not been specified.\n"
            << "Please correct the input deck.";
        }

        m_density_model = DensityModel::SpeciesMixture;
        m_h_species_densities.resize(m_nspecies, 0.);

      } else if (amrex::toLower(density_model) != "none") {

        reporter::Log(reporter::Error, __FILE__, __LINE__)
          << "Invalid solids density model specified.\n"
          << "Valid entries are: None; Mixture.\n"
          << "Please correct the input deck.";
      }

      for (int n(0); n < m_nspecies; n++) {
        const auto& names = species.names();
        auto it = std::find(names.begin(), names.end(), m_species_names[n]);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != names.end(),
                                         "Solids species missing in input");

        const auto pos = std::distance(names.begin(), it);

        m_species_IDs[n] = species.IDs(pos);
        m_MW_sn0[n] = species.MW_k(pos);

        if (m_density_model == DensityModel::SpeciesMixture) {

          m_h_species_densities[n] = species.density_k(pos);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_h_species_densities[n] > std::numeric_limits<amrex::Real>::min(),
                "Provided species density must be positive."
                );
        }
      }

    } else {

      if (pp_solids.contains("molecular_weight")) {
        m_MW_sn0.resize(1);
        pp_solids.get("molecular_weight", m_MW_sn0[0]);
      }
    }

    // Create the stoichiometric table for the fluid phase, that is associate
    // the total stoichiometric coefficient for each fluid species in each
    // reaction
    if (reactions.solve()) {

      const int lagrangian_nreactions = reactions.lagrangian_nreactions();

      // Allocate space for necessary data
      m_lagrangian_stoich_coeffs.resize(m_nspecies*lagrangian_nreactions, 0.);

      for (int n_s(0); n_s < m_nspecies; n_s++) {
        // Get the ID of the current species n_s
        const int species_id = m_species_IDs[n_s];

        // Loop over reactions to compute each contribution
        for (int q(0); q < lagrangian_nreactions; q++) {

          MFIXChemicalReaction* chem_reaction = reactions.get_lagrangian(q);
          const auto& chem_phases = chem_reaction->get_phases();

          const auto& reactants_IDs = chem_reaction->get_reactants_ids();
          const auto& reactants_chem_phases = chem_reaction->get_reactants_phases();
          const auto& reactants_coeffs = chem_reaction->get_reactants_coeffs();

          const auto& products_IDs = chem_reaction->get_products_ids();
          const auto& products_chem_phases = chem_reaction->get_products_phases();
          const auto& products_coeffs = chem_reaction->get_products_coeffs();

          // Do something only if reaction is lagrangian and contains
          // a solid compound
          if (chem_reaction->get_type() == ReactionType::Lagrangian &&
              std::find(chem_phases.begin(), chem_phases.end(), ChemicalPhase::Solid) != chem_phases.end()) {

            // Add reactant contribution (if any)
            {
              for (int pos(0); pos < reactants_IDs.size(); ++pos) {
                if (species_id == reactants_IDs[pos] && reactants_chem_phases[pos] == ChemicalPhase::Solid) {
                  m_lagrangian_stoich_coeffs[n_s*lagrangian_nreactions+q] += reactants_coeffs[pos];
                }
              }
            }

            // Add products contribution (if any)
            {
              for (int pos(0); pos < products_IDs.size(); ++pos) {
                if (species_id == products_IDs[pos] && products_chem_phases[pos] == ChemicalPhase::Solid) {
                  m_lagrangian_stoich_coeffs[n_s*lagrangian_nreactions+q] += products_coeffs[pos];
                }
              }
            }
          }
        }
      }

    }

    // Do not support multiple species for solids conductivities!
    // Always reads in the conductivity model from first solids phase
    std::string ks_model;
    if ( pp_solids.query("thermal_conductivity", ks_model) ) {

      if (amrex::toLower(ks_model).compare("constant") == 0) {

        m_thermal_conductivity_model = ThermalConductivityModel::Constant;

        m_kp_sn0.resize(1);
        pp_solids.query("thermal_conductivity.constant", m_kp_sn0[0]);

        if ( m_kp_sn0[0] <= 0. ) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Invalid DEM thermal conductivity: " << m_kp_sn0[0];
        }

        // Set flag to include conduction. The check on solving enthalpy
        // is probably overkill but what could it hurt?
        if (m_solve_enthalpy) { m_do_conduction = 1; }

        pp_solids.query("flpc", m_flpc);
        pp_solids.query("min_conduction_dist", m_rough);

      } else { // Error if the model is not constant.

        reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Invalid DEM thermal conductivity: " << ks_model;
      }
    }


    // TODO
    // move DEM initialization before Solids initialization and get this
    // information directly from DEM class
    int solve_dem = 1;
    {
      amrex::Vector<std::string> dem_names;
      amrex::ParmParse ppDEM("dem");
      ppDEM.queryarr("solve", dem_names);
      if (dem_names.size() > 1)
        if (amrex::toLower(dem_names[0]).compare("none") == 0 ||
            (dem_names[0]).compare("0") == 0)
          solve_dem = 0;
    }


    // Allocate m_h_parameters
    {
      m_h_parameters = new MFIXSolidsParms();

      m_h_parameters->m_nspecies = m_nspecies;
      m_h_parameters->m_solve_dem = solve_dem;

      // Species names
      if (m_h_species_names.size() > 0) {
        m_h_parameters->m_species_names = m_h_species_names.data();
        m_h_parameters->m_species_names_idxs = m_h_species_names_idxs.data();
        m_h_parameters->m_species_names_IDs = m_h_species_names_IDs.data();
      }

      // Molecular weights
      if (m_MW_sn0.size() > 0) {
        m_h_parameters->m_MW_sn = m_MW_sn0.data();
      }

      // Solid phase species densities
      if (m_h_species_densities.size() > 0) {
        m_h_parameters->m_species_densities = m_h_species_densities.data();
      }

      m_h_parameters->m_lagrangian_nreactions = reactions.lagrangian_nreactions();

      // Stoichiometric coefficients
      if (m_lagrangian_stoich_coeffs.size() > 0) {
        m_h_parameters->m_lagrangian_stoich_coeffs = m_lagrangian_stoich_coeffs.data();
      }

      //
      if (m_kp_sn0.size() > 0) {
        m_h_parameters->m_kp_sn = m_kp_sn0.data();
      }

      m_h_parameters->m_flpc = m_flpc;
      m_h_parameters->m_rough = m_rough;
      m_h_parameters->m_do_pfp_cond = m_do_conduction;

      m_h_parameters->m_thermal_conductivity_model = m_thermal_conductivity_model;
    }

    // Allocate m_d_parameters
#ifdef AMREX_USE_GPU
    {
      m_d_parameters = new MFIXSolidsParms();

      m_d_parameters->m_nspecies = m_nspecies;
      m_d_parameters->m_solve_dem = solve_dem;

      // Species names
      if (m_h_species_names.size() > 0) {
        m_d_species_names.resize(m_h_species_names.size());
        Gpu::copyAsync(Gpu::hostToDevice, m_h_species_names.begin(), m_h_species_names.end(), m_d_species_names.begin());
        m_d_parameters->m_species_names = m_d_species_names.data();

        m_d_species_names_idxs.resize(m_h_species_names_idxs.size());
        Gpu::copyAsync(Gpu::hostToDevice, m_h_species_names_idxs.begin(), m_h_species_names_idxs.end(), m_d_species_names_idxs.begin());
        m_d_parameters->m_species_names_idxs = m_d_species_names_idxs.data();

        m_d_species_names_IDs.resize(m_h_species_names_IDs.size());
        Gpu::copyAsync(Gpu::hostToDevice, m_h_species_names_IDs.begin(), m_h_species_names_IDs.end(), m_d_species_names_IDs.begin());
        m_d_parameters->m_species_names_IDs = m_d_species_names_IDs.data();
      }

      // Molecular weights
      if (m_MW_sn0.size() > 0) {
        m_d_MW_sn0.resize(m_MW_sn0.size());
        Gpu::copyAsync(Gpu::hostToDevice, m_MW_sn0.begin(), m_MW_sn0.end(), m_d_MW_sn0.begin());
        m_d_parameters->m_MW_sn = m_d_MW_sn0.data();
      }

      // Solid phase species densities
      if (m_h_species_densities.size() > 0) {
        m_d_species_densities.resize(m_h_species_densities.size());
        Gpu::copyAsync(Gpu::hostToDevice, m_h_species_densities.begin(), m_h_species_densities.end(), m_d_species_densities.begin());
        m_d_parameters->m_species_densities = m_d_species_densities.data();
      }

      m_d_parameters->m_lagrangian_nreactions = reactions.lagrangian_nreactions();

      // Stoichiometric coefficients
      if (m_lagrangian_stoich_coeffs.size() > 0) {
        m_d_lagrangian_stoich_coeffs.resize(m_lagrangian_stoich_coeffs.size());

        Gpu::copyAsync(Gpu::hostToDevice, m_lagrangian_stoich_coeffs.begin(),
            m_lagrangian_stoich_coeffs.end(), m_d_lagrangian_stoich_coeffs.begin());

        m_d_parameters->m_lagrangian_stoich_coeffs = m_d_lagrangian_stoich_coeffs.data();
      }

      //
      if (m_kp_sn0.size() > 0) {
        m_d_kp_sn0.resize(m_kp_sn0.size());
        Gpu::copyAsync(Gpu::hostToDevice, m_kp_sn0.begin(), m_kp_sn0.end(), m_d_kp_sn0.begin());
        m_d_parameters->m_kp_sn = m_d_kp_sn0.data();
      }

      m_d_parameters->m_flpc = m_flpc;
      m_d_parameters->m_rough = m_rough;
      m_d_parameters->m_do_pfp_cond = m_do_conduction;

      m_d_parameters->m_thermal_conductivity_model = m_thermal_conductivity_model;
    }
#endif

  }
}




long
SOLIDS_t::
calc_particles_from_vol (amrex::Real const a_vol,
                         amrex::Real const a_multi_particle_volume,
                         int const include_volfrac) const
{
  Real vfrac = (include_volfrac ? volfrac : 1.);

  Real np_r(0.);

  if (a_multi_particle_volume > 0.) {

    np_r = (vfrac * a_vol) / a_multi_particle_volume;

  } else if (diameter.is_constant()) {

    Real const dp = diameter.get_mean();
    Real const vol_p = (M_PI/6.0)*dp*dp*dp;

    np_r = (vfrac * a_vol) / vol_p;

  } else { //if (diameter.is_normal() || diameter.is_uniform() ) {

    //diameter.report_distribution(1.0e6);

    for ( int bin(0); bin < diameter.get_bins(); ++bin) {

      auto const dist_info = diameter.get_Vw_probability(bin);

      Real const dp = std::get<0>(dist_info);
      Real const vol_p = (M_PI/6.0)*dp*dp*dp;

      Real const probability = std::get<1>(dist_info);
      Real const vol_s = vfrac * a_vol * probability;

#if 0
        PrintToFile("dist", Print::AllProcs)
          << std::setw(3) << bin << "  "
          << std::scientific << std::setw(16) << std::setprecision(8) << dp << "   "
          << std::scientific << std::setw(16) << std::setprecision(8) << vol_s << "   "
          << std::scientific << std::setw(16) << std::setprecision(8) << vol_p << "   "
          << std::fixed      << std::setw( 8) << std::setprecision(4) << probability << "   "
          << std::fixed      << std::setw( 8) << std::setprecision(4) << (vol_s/vol_p) << "   "
          //<< std::setw(8) << static_cast<long>(vol_s / vol_p)
          << '\n';
#endif
      np_r += vol_s / vol_p;
    }

  }

  return static_cast<long>(np_r);
}
