#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_EBMultiFabUtil.H>

#include <mfix_reporter.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_fix_inputs.H>
#include <mfix_solids.H>

using namespace amrex;


AMREX_GPU_HOST_DEVICE
MFIXFluidParms::MFIXFluidParms ()
  : m_k_g(0.)
  , m_nspecies(0)
  , m_species_names(nullptr)
  , m_species_names_idxs(nullptr)
  , m_species_names_IDs(nullptr)
  , m_MW_gk(nullptr)
  , m_D_g(0)
  , m_lagrangian_nreactions(0)
  , m_eulerian_nreactions(0)
  , m_lagrangian_stoich_coeffs(nullptr)
  , m_eulerian_stoich_coeffs(nullptr)
{}


MFIXFluidPhase::MFIXFluidPhase ()
  : m_ntypes(0)
  , m_names(0)
  , m_solve(0)
  , m_solve_density(0)
  , m_solve_tracer(0)
  , m_ntracer(0)
  , m_trac_0(0)
  , m_solve_enthalpy(0)
  , m_k_g0(0)
  , m_thermodynamic_pressure(-1.0)
  , m_p_therm_defined(0)
  , m_solve_species(0)
  , m_species_names(0)
  , m_h_species_names(0)
  , m_h_species_names_idxs(0)
  , m_d_species_names(0)
  , m_d_species_names_idxs(0)
  , m_species_IDs(0)
  , m_nspecies(0)
  , m_MW_gk0(0)
  , m_d_MW_gk0(0)
  , m_D_g0(0)
  , m_is_a_mixture(0)
  , m_lagrangian_stoich_coeffs(0)
  , m_d_lagrangian_stoich_coeffs(0)
  , m_eulerian_stoich_coeffs(0)
  , m_d_eulerian_stoich_coeffs(0)
  , m_h_parameters(nullptr)
  , m_d_parameters(nullptr)
  , m_is_initialized(0)
{}

MFIXFluidPhase::~MFIXFluidPhase()
{
  if (m_h_parameters != nullptr)
    delete m_h_parameters;

  if (m_d_parameters != nullptr)
    delete m_d_parameters;
}

int MFIXFluidPhase::
Initialize (const MFIXSpecies& species,
            const MFIXReactions& reactions)
{
  ParmParse pp;

  pp.queryarr("fluid.solve", m_names);
  m_ntypes = m_names.size();

  // Abort if more than one fluid type because it's not yet implemented
  if (m_ntypes > 1) {
    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Multi-fluid models are not currently supported.\n"
      << "Please correct the input deck.";
    return 1;
  }

  // Disable the fluid solver if the fluid is defined as "" or "None"
  if ((m_ntypes > 0) && (toLower(m_names[0]) != "none")) {

    if (!species.isInitialized()) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "The species object passed to the fluid is uninitialized!";
      return 1;
    }

    if (!reactions.isInitialized()) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "The reactions object passed to the fluid is uninitialized!";
      return 1;
    }

    ParmParse ppFluid(m_names[0]);

    // Set the fluid solving flag to true
    m_solve = 1;

    // Get mfix global inputs ------------------------------------------------//
    pp.query("mfix.advect_density",  m_solve_density);
    pp.query("mfix.advect_enthalpy", m_solve_enthalpy);
    pp.query("mfix.solve_species",   m_solve_species);
    pp.query("mfix.advect_tracer",   m_solve_tracer);


    { // Constraint type

      // Fix old way for specifying drag model
      FixInputs fix("Oct. 2024");
      fix.swap<std::string>("mfix.constraint_type", "mfix.constraint");

      std::string constraint_str = "IncompressibleFluid";
      if(!pp.queryAdd("mfix.constraint", constraint_str)) {

        reporter::Log(reporter::Status)
          << "Fluid constraint not specified. Setting default:\n"
          << "  mfix.constraint = " << constraint_str;
      }

      if ( constraint.set(constraint_str) ) {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Unknown fluid constraint type: " << constraint_str << "\n"
          << "Please correct the input deck.";
        return 1;
      }

      int include_depdt(0);
      if (pp.query("mfix.constraint.include_depdt", include_depdt))
      { constraint.set_depdt(include_depdt); }

      if ( constraint.include_depdt() ) {
#if 1
        reporter::Log(reporter::Warning) << "Including d(epf)/dt in constraint.\n"
          << "This is option is under development. Consider yourself warned.";
#else
        reporter::Log(reporter::Error,__FILE__,__LINE__)
          << "Option to including d(epg)/dt in constraint is disabled.\n"
          << "This is option is under development and not in a working state.";
#endif
      }
    }

    // Fix old way for specifying molecular viscosity
    { FixInputs fix(1);

      fix.swap<std::string>("fluid.viscosity", "fluid.viscosity.molecular", 1);
      fix.swap<amrex::Real>("fluid.viscosity.constant", "fluid.viscosity.molecular.constant", 1);
    }

    // Get eddy viscosity inputs ----------------------------------//
    {
      std::string model("None");
      ppFluid.query("viscosity.eddy", model);
      std::string mlower(toLower(model));

      if (mlower == "smagorinsky-lilly") {
        m_eddy_visc_model = EddyViscModel::Smagorinsky;

        if (!ppFluid.query("viscosity.eddy.Smagorinsky-Lilly.constant", m_smagorinsky_constant) ) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Smagorinsky-Lilly constant is not defined!\n"
            << "Please correct the input deck.";
          return 1;
        }

      } else if (mlower == "wale") {
        m_eddy_visc_model = EddyViscModel::WALE;
        m_wale_constant = 0.325;

        if (!ppFluid.query("viscosity.eddy.WALE.constant", m_wale_constant) ) {
          reporter::Log(reporter::Info,__FILE__, __LINE__)
            << "WALE constant is not specified in the input deck\n"
            << "Using default " << m_wale_constant << "\n";
        }

      } else if (mlower == "usr") {
        m_eddy_visc_model = EddyViscModel::Usr;

      } else if (mlower == "none") {
        m_eddy_visc_model = EddyViscModel::None;

      } else {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Unknown eddy viscosity model: " << model << '\n'
          << "Please correct the input deck.";
        return 1;
      }
    } // end eddy viscosity inputs

    // Get suspension viscosity inputs ----------------------------------//
    {
      std::string model("None");
      ppFluid.query("viscosity.suspension", model);

      std::string mlower(toLower(model));
      if (mlower == "einstein") {

        m_susp_visc_model = SuspViscModel::Einstein;

      } else if (mlower == "brinkman") {

        m_susp_visc_model = SuspViscModel::Brinkman;

        if (!ppFluid.query("viscosity.suspension.Brinkman.constant",
                  m_brinkman_constant) ) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Brinkman viscosity constant is not defined!\n"
            << "Please correct the input deck.";
          return 1;
        }

      } else if (mlower == "roscoe") {

        m_susp_visc_model = SuspViscModel::Roscoe;

        if (!ppFluid.query("viscosity.suspension.Roscoe.c1", m_roscoe_c1) ) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Roscoe viscosity constant 1 is not defined!\n"
            << "Please correct the input deck.";
          return 1;
        }

        if (!ppFluid.query("viscosity.suspension.Roscoe.c2", m_roscoe_c2) ) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Roscoe viscosity constant 2 is not defined!\n"
            << "Please correct the input deck.";
          return 1;
        }

      } else if (mlower == "chenglaw") {

        m_susp_visc_model = SuspViscModel::ChengLaw;

        if (!ppFluid.query("viscosity.suspension.ChengLaw.constant", m_chenglaw_constant) ) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "ChengLaw viscosity constant is not defined!\n"
            << "Please correct the input deck.";
          return 1;
        }

      } else if (mlower == "sato") {

        m_susp_visc_model = SuspViscModel::Sato;
        m_sato_constant = 0.65;

        if (!ppFluid.query("viscosity.suspension.Sato.constant", m_sato_constant) ) {
          reporter::Log(reporter::Info,__FILE__, __LINE__)
            << "Sato constant is not specified in the input deck\n"
            << "Using default " << m_sato_constant;
        }

      } else if (mlower == "subramaniam") {

        m_susp_visc_model = SuspViscModel::Subramaniam;
        m_subramaniam_constant = 0.09;

        if (!ppFluid.query("viscosity.suspension.Subramaniam.constant", m_subramaniam_constant) ) {
          reporter::Log(reporter::Info,__FILE__, __LINE__)
            << "Subramaniam constant is not specified in the input deck\n"
            << "Using default " << m_subramaniam_constant;
        }

      } else if (mlower == "none") {

        m_susp_visc_model = SuspViscModel::None;

      } else {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Unknown suspension viscosity model: " << model << '\n'
          << "Please correct the input deck.";
        return 1;
      }

    } // end suspension viscosity inputs

    // Input for max effective viscosity factor
    {
      m_max_effective_viscosity_factor = 1.e6;

      if (!ppFluid.query("viscosity.max_effective_factor", m_max_effective_viscosity_factor) ) {
        reporter::Log(reporter::Info,__FILE__, __LINE__)
          << "Max effective viscosity factor not specified in the input deck\n"
          << "Using default " << m_max_effective_viscosity_factor;
      }
    } // end max effective viscosity factor input


    // Get specific heat and enthalpy inputs ---------------------------//
    if (m_solve_enthalpy) {
      if (props.specificHeat.define(m_names[0]) != 0)
      {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
           << "Failed to initialize specific heat!";
        return 1;
      }
    }

    // Molecular weight --------------------------------------------//
    if( !m_solve_species ) {

      int MW_reqd(0);

      if (constraint.isIdealGas()) { MW_reqd = 1; }

      if ( m_solve_enthalpy && (props.specificHeat.m_model == SpecificHeatModel::NASA7Polynomials)) {
          MW_reqd = 1;
      }

      if ( MW_reqd ) {

        m_MW_gk0.resize(1);

        if( !ppFluid.query("molecular_weight", m_MW_gk0[0]) ) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Fluid molecular weight was not specified!\n"
            << "Please correct the input deck.";
          return 1;
        }
        if ( m_MW_gk0[0] <= 0) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Invalid fluid molecular weight specified!\n"
            << "Please correct the input deck.";
          return 1;
        }
      }
    }

    if (m_solve_tracer) {

      m_ntracer = 1;
      pp.query("tracer", m_ntracer);

      if ( m_ntracer != 1) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Current tracer implementation only supports one tracer!";
      }

      // Scalar diffusive coefficient
      m_tracer_diff_coeff.resize(m_ntracer, 0.);
      pp.queryarr("tracer.diff_coeff", m_tracer_diff_coeff, 0, m_ntracer);

      Print() << "Scalar diffusion coefficients\n";
      for (int i(0); i < m_ntracer; i++)
      { Print() << "Tracer" << i << ": " << m_tracer_diff_coeff[i] << '\n'; }

      // Query fluid tracer initial value
      ppFluid.query("trac0", m_trac_0);
    } // solve tracer

    if (m_solve_species) {
      // Query fluid species
      ppFluid.queryarr("species", m_species_names);
    }

    if ((!m_solve_species) || (m_species_names.size() == 0)) {
      m_species_names.clear();
      m_h_species_names.clear();
      m_h_species_names_idxs.clear();
      m_solve_species = 0;
      m_nspecies = 0;
    } else if (toLower(m_species_names[0]) == "none") {
      m_species_names.clear();
      m_h_species_names.clear();
      m_h_species_names_idxs.clear();
      m_solve_species = 0;
      m_nspecies = 0;
    } else {
      m_solve_species = 1;
      m_nspecies = m_species_names.size();

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_species_names.size() > 0,
                                       "No input provided for fluid.names_names");

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nspecies <= species.nspecies(),
          "Fluid species number is higher than total species number");

      if (m_nspecies > 0) {

        std::map<std::string, int> ordered_species_names;
        int length(0);
        for (int n(0); n < m_nspecies; ++n) {
          const std::string name = m_species_names[n];
          ordered_species_names.insert(std::make_pair(name, n));
          length += name.size();
        }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nspecies == Long(ordered_species_names.size()),
            "There are fluid species defined multiple times. Check inputs");

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

    // Flag to determine if we want to solve the fluid as a mixture
    m_is_a_mixture = static_cast<int>(m_nspecies > 1);

    if (m_solve_enthalpy) {

      // Get thermal conductivity inputs -----------------------------//
      {
        std::string model;
        if (!ppFluid.query("thermal_conductivity", model)) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "The fluid thermal conductivity model was not specified!\n"
            << "Please correct the input deck.";
          return 1;
        }

        if (toLower(model) ==  "constant") {
          m_thermal_conductivity_model = ThermalConductivityModel::Constant;
          if (!ppFluid.query("thermal_conductivity.constant", m_k_g0)) {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Constant fluid thermal conductivity was not specified!\n"
              << "Please correct the input deck.";
            return 1;
          }

        } else {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Unknown fluid thermal conductivity model: " << model << '\n'
            << "Please correct the input deck.";
          return 1;
        }
      }

    }

    // Get fluid species m_parameters from species class
    if (m_solve_species) {

      m_species_IDs.resize(m_nspecies);
      m_MW_gk0.resize(m_nspecies);
      m_D_g0 = species.diffusivity();

      for (int n(0); n < m_nspecies; n++) {
        const auto& names = species.names();
        auto it = std::find(names.begin(), names.end(), m_species_names[n]);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != names.end(),
                                         "Fluid species missing in input");

        const auto pos = std::distance(names.begin(), it);

        m_species_IDs[n] = species.IDs(pos);
        m_MW_gk0[n] = species.MW_k(pos);
      }

    }

    // Create the stoichiometric table for the fluid phase, that is associate
    // the total stoichiometric coefficient for each fluid species in each
    // reaction
    if (reactions.solve()) {

      const int lagrangian_nreactions = reactions.lagrangian_nreactions();
      const int eulerian_nreactions = reactions.eulerian_nreactions();

      // Allocate space for necessary data
      m_lagrangian_stoich_coeffs.resize(m_nspecies*lagrangian_nreactions, 0.);

      for (int n_g(0); n_g < m_nspecies; n_g++) {
        // Get the ID of the current species n_g
        const int species_id = m_species_IDs[n_g];

        // Loop over reactions to compute each contribution
        for (int q(0); q < lagrangian_nreactions; q++) {

          MFIXChemicalReaction* chem_reaction = reactions.get_lagrangian(q);
          const auto& phases = chem_reaction->get_phases();

          const auto& reactants_IDs = chem_reaction->get_reactants_ids();
          const auto& reactants_phases = chem_reaction->get_reactants_phases();
          const auto& reactants_coeffs = chem_reaction->get_reactants_coeffs();

          const auto& products_IDs = chem_reaction->get_products_ids();
          const auto& products_phases = chem_reaction->get_products_phases();
          const auto& products_coeffs = chem_reaction->get_products_coeffs();

          // Do something only if reaction is lagrangian and contains a fluid
          // compound
          if (chem_reaction->get_type() == ReactionType::Lagrangian &&
              std::find(phases.begin(), phases.end(), ChemicalPhase::Fluid) != phases.end()) {

            // Add reactant contribution (if any)
            for (int pos(0); pos < reactants_IDs.size(); ++pos) {
              if (species_id == reactants_IDs[pos] &&
                  reactants_phases[pos] == ChemicalPhase::Fluid) {
                m_lagrangian_stoich_coeffs[n_g*lagrangian_nreactions+q] += reactants_coeffs[pos];
              }
            }

            // Add products contribution (if any)
            for (int pos(0); pos < products_IDs.size(); ++pos) {
              if (species_id == products_IDs[pos] &&
                  products_phases[pos] == ChemicalPhase::Fluid) {
                m_lagrangian_stoich_coeffs[n_g*lagrangian_nreactions+q] += products_coeffs[pos];
              }
            }
          }
        }
      }


      // Allocate space for necessary data
      m_eulerian_stoich_coeffs.resize(m_nspecies*eulerian_nreactions, 0.);

      for (int n_g(0); n_g < m_nspecies; n_g++) {
        // Get the ID of the current species n_g
        const int species_id = m_species_IDs[n_g];

        // Loop over reactions to compute each contribution
        for (int q(0); q < eulerian_nreactions; q++) {

          MFIXChemicalReaction* chem_reaction = reactions.get_eulerian(q);
          const auto& phases = chem_reaction->get_phases();

          const auto& reactants_IDs = chem_reaction->get_reactants_ids();
          const auto& reactants_phases = chem_reaction->get_reactants_phases();
          const auto& reactants_coeffs = chem_reaction->get_reactants_coeffs();

          const auto& products_IDs = chem_reaction->get_products_ids();
          const auto& products_phases = chem_reaction->get_products_phases();
          const auto& products_coeffs = chem_reaction->get_products_coeffs();

          // Do something only if reaction is eulerian and contains a fluid
          // compound
          if (chem_reaction->get_type() == ReactionType::Eulerian &&
              std::find(phases.begin(), phases.end(), ChemicalPhase::Fluid) != phases.end()) {

            // Add reactant contribution (if any)
            for (int pos(0); pos < reactants_IDs.size(); ++pos) {
              if (species_id == reactants_IDs[pos] &&
                  reactants_phases[pos] == ChemicalPhase::Fluid) {
                m_eulerian_stoich_coeffs[n_g*eulerian_nreactions+q] += reactants_coeffs[pos];
              }
            }

            // Add products contribution (if any)
            for (int pos(0); pos < products_IDs.size(); ++pos) {
              if (species_id == products_IDs[pos] &&
                  products_phases[pos] == ChemicalPhase::Fluid) {
                m_eulerian_stoich_coeffs[n_g*eulerian_nreactions+q] += products_coeffs[pos];
              }
            }
          }
        }
      }

    }


    m_p_therm_defined = ppFluid.query("thermodynamic_pressure", m_thermodynamic_pressure);

    if (constraint.isIncompressibleFluid() && m_p_therm_defined) {

      amrex::Warning("When the incompressible fluid constraint is selected, "
          "the fluid thermodynamic pressure input will be ignored");
    }

    if (constraint.isIdealGasClosedSystem() && (!m_p_therm_defined)) {

      amrex::Abort("When the idealgas closedsystem constraint is selected, "
          "the fluid thermodynamic pressure input must be provided");
    }

    // Check on inputs in case of Ideal Gas EOS
    if ( constraint.isIdealGas() ) {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_MW_gk0.size() > 0, "Inputs error: fluid molecular_weight not provided");

      for (size_t i(0); i < m_MW_gk0.size(); ++i) {

        if (m_MW_gk0[i] < 1.e-15) {

          Print() << "Invalid molecular weight for species " << m_species_names[i] << "\n";
          amrex::Abort("Inputs error");
        }
      }
    }

    // Allocate m_h_parameters
    {
      m_h_parameters = new MFIXFluidParms();

      m_h_parameters->m_k_g = m_k_g0;
      m_h_parameters->m_nspecies = m_nspecies;

      // Species names
      if (m_h_species_names.size() > 0) {
        m_h_parameters->m_species_names = m_h_species_names.data();
        m_h_parameters->m_species_names_idxs = m_h_species_names_idxs.data();
        m_h_parameters->m_species_names_IDs = m_h_species_names_IDs.data();
      }

      if (m_MW_gk0.size() > 0) {
        m_h_parameters->m_MW_gk = m_MW_gk0.data();
      }

      m_h_parameters->m_D_g = m_D_g0;

      m_h_parameters->m_lagrangian_nreactions = reactions.lagrangian_nreactions();
      m_h_parameters->m_eulerian_nreactions = reactions.eulerian_nreactions();

      if (m_lagrangian_stoich_coeffs.size() > 0) {
        m_h_parameters->m_lagrangian_stoich_coeffs = m_lagrangian_stoich_coeffs.data();
      }

      if (m_eulerian_stoich_coeffs.size() > 0) {
        m_h_parameters->m_eulerian_stoich_coeffs = m_eulerian_stoich_coeffs.data();
      }

    }

    // Allocate m_d_parameters
#ifdef AMREX_USE_GPU
    {
      m_d_parameters = new MFIXFluidParms();

      m_d_parameters->m_k_g = m_k_g0;
      m_d_parameters->m_nspecies = m_nspecies;

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

      if (m_MW_gk0.size() > 0) {
        m_d_MW_gk0.resize(m_MW_gk0.size());
        Gpu::copyAsync(Gpu::hostToDevice, m_MW_gk0.begin(), m_MW_gk0.end(), m_d_MW_gk0.begin());
        m_d_parameters->m_MW_gk = m_d_MW_gk0.data();
      }

      m_d_parameters->m_D_g = m_D_g0;

      m_d_parameters->m_lagrangian_nreactions = reactions.lagrangian_nreactions();
      m_d_parameters->m_eulerian_nreactions = reactions.eulerian_nreactions();

      if (m_lagrangian_stoich_coeffs.size() > 0) {
        m_d_lagrangian_stoich_coeffs.resize(m_lagrangian_stoich_coeffs.size());

        Gpu::copyAsync(Gpu::hostToDevice, m_lagrangian_stoich_coeffs.begin(),
            m_lagrangian_stoich_coeffs.end(), m_d_lagrangian_stoich_coeffs.begin());

        m_d_parameters->m_lagrangian_stoich_coeffs = m_d_lagrangian_stoich_coeffs.data();
      }

      if (m_eulerian_stoich_coeffs.size() > 0) {
        m_d_eulerian_stoich_coeffs.resize(m_eulerian_stoich_coeffs.size());

        Gpu::copyAsync(Gpu::hostToDevice, m_eulerian_stoich_coeffs.begin(),
            m_eulerian_stoich_coeffs.end(), m_d_eulerian_stoich_coeffs.begin());

        m_d_parameters->m_eulerian_stoich_coeffs = m_d_eulerian_stoich_coeffs.data();
      }

    }
#endif

    // Get molecular viscosity model inputs ---------------------------//
    {
      if (props.molWeight.define(m_names[0]) != 0)
      {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
           << "Failed to initialize molecular weight!";
        return 1;
      }

      if (props.molViscosity.define(m_names[0]) != 0)
      {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
           << "Failed to initialize molecular viscosity!";
        return 1;
      }
    }

  }

  // Flag for initialization
  m_is_initialized = 1;

  return 0;
}

Real
FLUID_t::get_density (amrex::Real a_time) const
{

  if ( density.is_defined()) {

    Vector<Real> v = density(a_time);
    return v.at(0);

  } else { // density not defined, compute from ideal gas

    if (fluid->constraint.isIdealGas()) {

      const auto& fluid_parms = fluid->parameters<RunOn::Host>();

      amrex::Real MW_g(0.);

      if (fluid->isMixture()) {
        for (int n(0); n < fluid->nspecies(); n++) {
          MW_g += this->get_species(n, a_time) / fluid_parms.get_MW_gk(n);
        }

        MW_g = 1./MW_g;
      } else {
        MW_g = fluid_parms.get_MW_g();
      }

      return (fluid->thermodynamic_pressure() * MW_g) /
             (MFIXFluidPhase::R * this->get_temperature(a_time));

    } else {

      amrex::Abort("Error. How did we even arrive here?");

    }
  }

  return std::numeric_limits<Real>::max();
}


Real
FLUID_t::
get_species (int a_n, Real a_time) const
{
  AMREX_ASSERT(a_n < species.size());
  if ( !species[a_n].is_defined()) { return 0.; }
  Vector<Real> v = species[a_n](a_time);
  return v.at(0);
}


RealVect
FLUID_t::get_velocity (amrex::Real a_time) const
{
  if ( !velocity.is_defined()) {
    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Requested velocity magnitude not defined.";
  }

  Vector<Real> v = velocity(a_time);
  if (v.size() == 3) { return RealVect(v[0], v[1], v[2]); }
  else {
    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Requested BC velocity has invalid size.";
  }
  return RealVect(std::numeric_limits<Real>::max());
}


Real
FLUID_t::get_velocity_mag (amrex::Real a_time) const
{
  if ( velocity.is_defined()) {
    Vector<Real> v = velocity(a_time);
    if (v.size() == 1) {
      return v[0];
    } else {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Requested BC velocity mag has invalid size!";
    }
  } else {
    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Requested velocity magnitude not defined.";
  }
  return std::numeric_limits<Real>::max();
}


Real
FLUID_t::get_volflow (amrex::Real a_time) const
{
  if ( volflow.is_defined()) {
    Vector<Real> v = volflow(a_time);
    return v.at(0);
  } else {
    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Requested volflow not defined.";
  }
  return std::numeric_limits<Real>::max();
}


Real
FLUID_t::get_massflow (amrex::Real a_time) const
{
  if ( massflow.is_defined()) {
    Vector<Real> v = massflow(a_time);
    return v.at(0);
  } else {
    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Requested massflow not defined.";
  }
  return std::numeric_limits<Real>::max();
}


Real
FLUID_t::get_temperature (amrex::Real a_time) const
{
  if ( temperature.is_defined()) {
    Vector<Real> v = temperature(a_time);
    return v.at(0); }
  else {
    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Requested temperature not defined.";
  }

  return std::numeric_limits<Real>::max();
}
