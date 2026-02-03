#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>

#include <mfix_reactions.H>
#include <mfix_reporter.H>

#include <string>
#include <sstream>
#include <cctype>
#include <iterator>
#include <regex>
#include <algorithm>

using namespace amrex;


namespace chemistry_aux {

std::string ltrim(const std::string& s) {
  return std::regex_replace(s, std::regex("^\\s+"), std::string(""));
}


std::string rtrim(const std::string& s) {
  return std::regex_replace(s, std::regex("\\s+$"), std::string(""));
}


std::string trim(const std::string& s) {
  return ltrim(rtrim(s));
}

} // end namespace chemistry_aux


MFIXReactions::MFIXReactions()
  : m_solve(0)
  , m_lagrangian_nreactions(0)
  , m_eulerian_nreactions(0)
  , m_reactions_names(0)
  , m_reactions_equations(0)
  , m_h_eulerian_reactions_names(0)
  , m_h_eulerian_reactions_names_idxs(0)
  , m_h_eulerian_reactions_names_IDs(0)
  , m_d_eulerian_reactions_names(0)
  , m_d_eulerian_reactions_names_idxs(0)
  , m_d_eulerian_reactions_names_IDs(0)
  , m_h_lagrangian_reactions_names(0)
  , m_h_lagrangian_reactions_names_idxs(0)
  , m_h_lagrangian_reactions_names_IDs(0)
  , m_d_lagrangian_reactions_names(0)
  , m_d_lagrangian_reactions_names_idxs(0)
  , m_d_lagrangian_reactions_names_IDs(0)
  , m_reactions_IDs(0)
  , m_nreactions(0)
  , m_is_initialized(0)
  , m_h_eulerian_parameters(nullptr)
  , m_d_eulerian_parameters(nullptr)
  , m_h_lagrangian_parameters(nullptr)
  , m_d_lagrangian_parameters(nullptr)
  , m_eulerian_chemical_reactions(0)
  , m_lagrangian_chemical_reactions(0)
{}


MFIXReactions::~MFIXReactions()
{
  if (m_solve) {

    for (int n(0); n < m_lagrangian_nreactions; n++)
      delete m_lagrangian_chemical_reactions[n];

    for (int n(0); n < m_eulerian_nreactions; n++)
      delete m_eulerian_chemical_reactions[n];
  }

  if (m_h_eulerian_parameters != nullptr)
    delete m_h_eulerian_parameters;

  if (m_d_eulerian_parameters != nullptr)
    delete m_d_eulerian_parameters;

  if (m_h_lagrangian_parameters != nullptr)
    delete m_h_lagrangian_parameters;

  if (m_d_lagrangian_parameters != nullptr)
    delete m_d_lagrangian_parameters;
}


// Initialization: read input parameters and set up reactions
void MFIXReactions::Initialize (const MFIXSpecies& species)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.isInitialized(),
      "Species not initialized. Can't initialize reactions before species initialization");

  // Set the initialization flag
  m_is_initialized = 1;

  amrex::ParmParse pp("chemistry");

  if (pp.contains("solve")) {

    // Get the list of reactions names to define chem equations
    pp.getarr("solve", m_reactions_names);

    m_nreactions = m_reactions_names.size();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nreactions > 0, "No input provided for chemistry.solve");

    // Disable the species solver if the species are defined as "None" (case
    // insensitive) or 0
    if (amrex::toLower(m_reactions_names[0]).compare("none") == 0 ||
        (m_reactions_names[0]).compare("0") == 0)
    {

      m_solve = 0;
      m_nreactions = 0;

    } else {

      {
        std::set<std::string> ordered_reactions_names;

        for (auto elem: m_reactions_names)
          ordered_reactions_names.insert(elem);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nreactions == Long(ordered_reactions_names.size()),
            "There are reactions defined multiple times. Check inputs");
      }

      m_solve = 1;

      m_reactions_equations.resize(m_nreactions);
      m_reactions_IDs.resize(m_nreactions);

      for (int n(0); n < m_nreactions; n++) {
        // Get the reation equation relative to given reaction name
        Vector<std::string> equation(0);
        pp.getarr((m_reactions_names[n]+".reaction"), equation);

        for (auto elem: equation)
          m_reactions_equations[n] += elem;

        m_reactions_IDs[n] = n;
      }


      std::map<std::string, int> ordered_eulerian_reactions_names;
      std::map<std::string, int> ordered_lagrangian_reactions_names;
      int eulerian_length(0);
      int lagrangian_length(0);

      // Loop over reactions
      for (int n(0); n < m_nreactions; ++n) {

        const std::string& name = m_reactions_names[n];
        const std::string& equation = m_reactions_equations[n];

        auto* reaction = new MFIXChemicalReaction(name, equation, species);

        if (reaction->get_type() == ReactionType::Eulerian) {

          m_eulerian_chemical_reactions.push_back(std::move(reaction));
          ordered_eulerian_reactions_names.insert(std::make_pair(name, m_eulerian_nreactions));
          eulerian_length += name.size();
          m_eulerian_nreactions++;

        } else if (reaction->get_type() == ReactionType::Lagrangian) {

          m_lagrangian_chemical_reactions.push_back(std::move(reaction));
          ordered_lagrangian_reactions_names.insert(std::make_pair(name, m_lagrangian_nreactions));
          lagrangian_length += name.size();
          m_lagrangian_nreactions++;

        } else {
          amrex::Abort("Error");
        }
      }


      // Set eulerian aux variables
      m_h_eulerian_reactions_names.clear();
      m_h_eulerian_reactions_names.reserve(eulerian_length+1);

      m_h_eulerian_reactions_names_idxs.clear();
      m_h_eulerian_reactions_names_idxs.reserve(m_eulerian_nreactions+1);

      m_h_eulerian_reactions_names_IDs.clear();
      m_h_eulerian_reactions_names_IDs.reserve(m_eulerian_nreactions);

      eulerian_length = 0;
      for (auto elem: ordered_eulerian_reactions_names) {
        const std::string& name = get<0>(elem);
        for (auto character: name) {
          m_h_eulerian_reactions_names.push_back(character);
        }
        m_h_eulerian_reactions_names_idxs.push_back(eulerian_length);
        m_h_eulerian_reactions_names_IDs.push_back(get<1>(elem));
        eulerian_length += name.size();
      }
      m_h_eulerian_reactions_names.push_back('\0');
      m_h_eulerian_reactions_names_idxs.push_back(eulerian_length);


      // Set lagrangian aux variables
      m_h_lagrangian_reactions_names.clear();
      m_h_lagrangian_reactions_names.reserve(lagrangian_length+1);

      m_h_lagrangian_reactions_names_idxs.clear();
      m_h_lagrangian_reactions_names_idxs.reserve(m_lagrangian_nreactions+1);

      m_h_lagrangian_reactions_names_IDs.clear();
      m_h_lagrangian_reactions_names_IDs.reserve(m_lagrangian_nreactions);

      lagrangian_length = 0;
      for (auto elem: ordered_lagrangian_reactions_names) {
        const std::string& name = get<0>(elem);
        for (auto character: name) {
          m_h_lagrangian_reactions_names.push_back(character);
        }
        m_h_lagrangian_reactions_names_idxs.push_back(lagrangian_length);
        m_h_lagrangian_reactions_names_IDs.push_back(get<1>(elem));
        lagrangian_length += name.size();
      }
      m_h_lagrangian_reactions_names.push_back('\0');
      m_h_lagrangian_reactions_names_idxs.push_back(lagrangian_length);

      // Flags for checking if a lagrangian reaction has both fluid and solids
      m_h_lagrangian_reactions_have_fluid_and_solids.clear();
      m_h_lagrangian_reactions_have_fluid_and_solids.reserve(m_lagrangian_nreactions);

      for (int i(0); i < m_lagrangian_nreactions; ++i) {
        auto const& reaction = m_lagrangian_chemical_reactions[i];
        m_h_lagrangian_reactions_have_fluid_and_solids.push_back(reaction->has_fluid_and_solids());
      }
    }
  }

  // Allocate m_h_parameters
  {
    m_h_eulerian_parameters = new MFIXReactionsParms(m_eulerian_nreactions);

    if (m_eulerian_nreactions > 0) {
      m_h_eulerian_parameters->m_reactions_names = m_h_eulerian_reactions_names.data();
      m_h_eulerian_parameters->m_reactions_names_idxs = m_h_eulerian_reactions_names_idxs.data();
      m_h_eulerian_parameters->m_reactions_names_IDs = m_h_eulerian_reactions_names_IDs.data();
      m_h_eulerian_parameters->m_reactions_have_fluid_and_solids = nullptr;
    }

    m_h_lagrangian_parameters = new MFIXReactionsParms(m_lagrangian_nreactions);

    if (m_lagrangian_nreactions > 0) {
      m_h_lagrangian_parameters->m_reactions_names = m_h_lagrangian_reactions_names.data();
      m_h_lagrangian_parameters->m_reactions_names_idxs = m_h_lagrangian_reactions_names_idxs.data();
      m_h_lagrangian_parameters->m_reactions_names_IDs = m_h_lagrangian_reactions_names_IDs.data();
      m_h_lagrangian_parameters->m_reactions_have_fluid_and_solids = m_h_lagrangian_reactions_have_fluid_and_solids.data();
    }
  }

  // Allocate m_d_parameters
#ifdef AMREX_USE_GPU
  {
    m_d_eulerian_parameters = new MFIXReactionsParms(m_eulerian_nreactions);

    if (m_eulerian_nreactions > 0) {
      m_d_eulerian_reactions_names.resize(m_h_eulerian_reactions_names.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_h_eulerian_reactions_names.begin(),
                                        m_h_eulerian_reactions_names.end(),
                                        m_d_eulerian_reactions_names.begin());
      m_d_eulerian_parameters->m_reactions_names = m_d_eulerian_reactions_names.data();

      m_d_lagrangian_reactions_names_idxs.resize(m_h_lagrangian_reactions_names_idxs.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_h_lagrangian_reactions_names_idxs.begin(),
                                        m_h_lagrangian_reactions_names_idxs.end(),
                                        m_d_lagrangian_reactions_names_idxs.begin());
      m_d_eulerian_parameters->m_reactions_names_idxs = m_d_lagrangian_reactions_names_idxs.data();

      m_d_lagrangian_reactions_names_IDs.resize(m_h_lagrangian_reactions_names_IDs.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_h_lagrangian_reactions_names_IDs.begin(),
                                        m_h_lagrangian_reactions_names_IDs.end(),
                                        m_d_lagrangian_reactions_names_IDs.begin());
      m_d_eulerian_parameters->m_reactions_names_IDs = m_d_lagrangian_reactions_names_IDs.data();
    }

    m_d_lagrangian_parameters = new MFIXReactionsParms(m_lagrangian_nreactions);

    if (m_lagrangian_nreactions > 0) {
      m_d_lagrangian_reactions_names.resize(m_h_lagrangian_reactions_names.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_h_lagrangian_reactions_names.begin(),
                                        m_h_lagrangian_reactions_names.end(),
                                        m_d_lagrangian_reactions_names.begin());
      m_d_lagrangian_parameters->m_reactions_names = m_d_lagrangian_reactions_names.data();

      m_d_lagrangian_reactions_names_idxs.resize(m_h_lagrangian_reactions_names_idxs.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_h_lagrangian_reactions_names_idxs.begin(),
                                        m_h_lagrangian_reactions_names_idxs.end(),
                                        m_d_lagrangian_reactions_names_idxs.begin());
      m_d_lagrangian_parameters->m_reactions_names_idxs = m_d_lagrangian_reactions_names_idxs.data();

      m_d_lagrangian_reactions_names_IDs.resize(m_h_lagrangian_reactions_names_IDs.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_h_lagrangian_reactions_names_IDs.begin(),
                                        m_h_lagrangian_reactions_names_IDs.end(),
                                        m_d_lagrangian_reactions_names_IDs.begin());
      m_d_lagrangian_parameters->m_reactions_names_IDs = m_d_lagrangian_reactions_names_IDs.data();

      m_d_lagrangian_reactions_have_fluid_and_solids.resize(m_h_lagrangian_reactions_have_fluid_and_solids.size());
      Gpu::copyAsync(Gpu::hostToDevice, m_h_lagrangian_reactions_have_fluid_and_solids.begin(),
                                        m_h_lagrangian_reactions_have_fluid_and_solids.end(),
                                        m_d_lagrangian_reactions_have_fluid_and_solids.begin());
      m_d_lagrangian_parameters->m_reactions_have_fluid_and_solids = m_d_lagrangian_reactions_have_fluid_and_solids.data();
    }
  }
#endif
}


// REACTION_T Class Constructor
MFIXChemicalReaction::MFIXChemicalReaction (const std::string& name,
                                            const std::string& reaction,
                                            const MFIXSpecies& species)
  : m_type(ReactionType::Undefined)
  , m_name(name)
  , m_formula(reaction)
  , m_phases(0)
  , m_reactants(0)
  , m_reactants_IDs(0)
  , m_reactants_coeffs(0)
  , m_reactants_phases(0)
  , m_products(0)
  , m_products_IDs(0)
  , m_products_coeffs(0)
  , m_products_phases(0)
  , m_mass_balance_tolerance(1.e-12)
{
  ParmParse pp_reactions("chemistry");
  pp_reactions.query("mass_balance_tolerance", m_mass_balance_tolerance);

  parse_reaction(species);
}


std::string
MFIXChemicalReaction::parse_reactants(const std::string& formula)
{
  {
    std::size_t pos = formula.find("-->");
    if(pos != std::string::npos)
      return chemistry_aux::trim(formula.substr(0,pos));
  }
  {
    std::size_t pos = formula.find("<--");
    if(pos != std::string::npos)
      return chemistry_aux::trim(formula.substr(pos+3, formula.size()-(pos+3)));
  }
  {
    std::size_t pos = formula.find("<=>");
    if(pos != std::string::npos) {
      //return trim(formula.substr(0,pos)) + "+" +
      //  trim(formula.substr(pos+3, formula.size()-(pos+3)));
      amrex::Abort("Not implemented");
      return std::string();
    }
  }

  amrex::Abort("Error: wrong format of chemical equation");
  return formula;
}


std::string
MFIXChemicalReaction::parse_products(const std::string& formula)
{
  {
    std::size_t pos = formula.find("<--");
    if(pos != std::string::npos)
      return chemistry_aux::trim(formula.substr(0,pos));
  }
  {
    std::size_t pos = formula.find("-->");
    if(pos != std::string::npos)
      return chemistry_aux::trim(formula.substr(pos+3, formula.size()-(pos+3)));
  }
  {
    std::size_t pos = formula.find("<=>");
    if(pos != std::string::npos) {
      //return trim(formula.substr(0,pos)) + "+" +
      //  trim(formula.substr(pos+3, formula.size()-(pos+3)));
      amrex::Abort("Not implemented");
      return std::string();
    }
  }

  amrex::Abort("Error: wrong format of chemical equation");
  return formula;
}


void
MFIXChemicalReaction::parse_stoichiometric_data(const std::string& s,
                                                amrex::Vector<std::string>& compounds,
                                                amrex::Vector<int>& compounds_id,
                                                amrex::Vector<amrex::Real>& coefficients,
                                                amrex::Vector<int>& phases,
                                                const MFIXSpecies& species)
{
  std::string formula(chemistry_aux::trim(s));
  std::replace(formula.begin(), formula.end(), '+', ' ');

  std::istringstream iss(formula);

  std::vector<std::string> stoichiometry((std::istream_iterator<std::string>(iss)),
                                          std::istream_iterator<std::string>());

  // Number of compounds
  const int nc = stoichiometry.size();

  compounds_id.clear();
  compounds.clear();
  coefficients.clear();
  phases.clear();

  compounds_id.resize(nc, -1);
  compounds.resize(nc, std::string());
  coefficients.resize(nc, 0.);
  phases.resize(nc, -1);

  for(int n(0); n < nc; n++)
  {
    std::string single_compound = stoichiometry[n];

    std::size_t pos(0);
    while (std::isdigit(single_compound.at(pos)) || single_compound.at(pos) == '.')
      pos++;

    if(pos == 0)
      coefficients[n] = 1.;
    else {
      coefficients[n] = std::stod(single_compound.substr(0, pos));
      single_compound.erase(0, pos);
    }

    pos = std::min(single_compound.find("(g)"), single_compound.find("(s)"));
    if (pos != std::string::npos) {
      std::string loc_phase = single_compound.substr(pos, pos+3);
      if (loc_phase.compare("(g)") == 0)
        phases[n] = ChemicalPhase::Fluid;
      else if (loc_phase.compare("(s)") == 0)
        phases[n] = ChemicalPhase::Solid;
      else
        amrex::Abort("Error: unrecognized phase in stoichiometric equation");

      single_compound.erase(pos, pos+3);
    }
    else {
      amrex::Abort("Error: wrong format for chemical equation");
    }

    compounds[n] = chemistry_aux::trim(single_compound);

    const auto& names = species.names();
    auto it = std::find(names.begin(), names.end(), compounds[n]);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != names.end(),
        "Error: species " + compounds[n] + " is present in a reaction but " +
        "it is not present in the list of declared species");

    const auto it_pos = std::distance(names.begin(), it);

    compounds_id[n] = species.IDs(it_pos);
  }
}


void
MFIXChemicalReaction::check_mass_balance (const MFIXSpecies& species)
{
  AMREX_ASSERT(m_reactants_IDs.size() == m_reactants_coeffs.size());
  AMREX_ASSERT(m_products_IDs.size() == m_products_coeffs.size());

  Real reactants_mass(0.);
  Real products_mass(0.);

  for (int n(0); n < m_reactants_IDs.size(); ++n) {
    int ID = m_reactants_IDs[n];
    reactants_mass += (-1.*m_reactants_coeffs[n])*species.MW_k(ID);
  }

  for (int n(0); n < m_products_IDs.size(); ++n) {
    int ID = m_products_IDs[n];
    products_mass += m_products_coeffs[n]*species.MW_k(ID);
  }

  Real diff_mass = std::abs(reactants_mass-products_mass);

  if(diff_mass > m_mass_balance_tolerance) {

    Print() << "\nUnbalanced reaction: " << m_formula << "\n\n";

    Print() << "Reactants mass balance: " << reactants_mass << "\n"
            << "Products mass balance: " << products_mass << "\n\n";

    Print() << "Mass balance difference: " << diff_mass << "\n"
            << "Mass balance tolerance: " << m_mass_balance_tolerance << "\n\n";

    amrex::Abort("Fix inputs either in species mass fractions or chemical reactions");
  }

  return;
}


void
MFIXChemicalReaction::parse_reaction(const MFIXSpecies& species)
{
  // remove spaces from reaction reaction_formula string
  std::string formula(m_formula);
  formula.erase(std::remove(formula.begin(), formula.end(), ' '), formula.end());

  // Clear reactants containers
  m_reactants.clear();
  m_reactants_coeffs.clear();
  m_reactants_phases.clear();

  // Clear products containers
  m_products.clear();
  m_products_coeffs.clear();
  m_products_phases.clear();

  // Get the reaction part of the formula
  std::string reaction_part = parse_reactants(formula);

  // Get the products part of the formula
  std::string production_part = parse_products(formula);

  // Get the reaction part stoichiometric coefficients, elements and phases
  parse_stoichiometric_data(reaction_part, m_reactants, m_reactants_IDs,
                            m_reactants_coeffs, m_reactants_phases, species);

  // Multiply reaction coefficients by -1
  for (int i(0); i < m_reactants_coeffs.size(); ++i)
    m_reactants_coeffs[i] *= -1;

  // Get the production part stoichiometric coefficients, elements and phases
  parse_stoichiometric_data(production_part, m_products, m_products_IDs,
                            m_products_coeffs, m_products_phases, species);

  check_mass_balance(species);

  // Fill m_phases with all the phases found in reactants
  for (const int phase: m_reactants_phases)
    if (std::find(m_phases.begin(), m_phases.end(), phase) == m_phases.end())
      m_phases.push_back(phase);

  // Add to m_phases all the phases found in products
  for (const int phase: m_products_phases)
    if (std::find(m_phases.begin(), m_phases.end(), phase) == m_phases.end())
      m_phases.push_back(phase);

  // Set the type of reaction, whether eulerian or lagrangian
  if (m_phases.size() == 1) {

    m_has_fluid_and_solids = 0;

    if (m_phases[0] == ChemicalPhase::Fluid) {
      m_type = ReactionType::Eulerian;
    } else if (m_phases[0] == ChemicalPhase::Solid) {
      m_type = ReactionType::Lagrangian;
    } else {
      amrex::Abort("Error: unrecognized phase");
    }

  } else if (m_phases.size() >= 2) {
    m_type = ReactionType::Lagrangian;

    int has_fluid(0);
    int has_solids(0);

    amrex::Vector<int> IDs = {};
    amrex::Vector<amrex::Real> coeffs = {};
    amrex::Vector<int> phases = {};

    IDs.insert(IDs.end(), m_reactants_IDs.begin(), m_reactants_IDs.end());
    IDs.insert(IDs.end(), m_products_IDs.begin(), m_products_IDs.end());

    coeffs.insert(coeffs.end(), m_reactants_coeffs.begin(), m_reactants_coeffs.end());
    coeffs.insert(coeffs.end(), m_products_coeffs.begin(), m_products_coeffs.end());

    phases.insert(phases.end(), m_reactants_phases.begin(), m_reactants_phases.end());
    phases.insert(phases.end(), m_products_phases.begin(), m_products_phases.end());

    AMREX_ASSERT(IDs.size() == coeffs.size());
    AMREX_ASSERT(IDs.size() == phases.size());

    for (int i(0); i < IDs.size(); ++i) {
      int ID_i = IDs[i];
      amrex::Real coeff_i = coeffs[i];
      int phase_i = phases[i];

      for (int j(i+1); j < IDs.size(); ++j) {
        int ID_j = IDs[j];
        amrex::Real coeff_j = coeffs[j];
        int phase_j = phases[j];

        if ((ID_i == ID_j) && (phase_i == phase_j)) {
          coeff_i += coeff_j;
        }
      }

      if (!amrex::almostEqual(coeff_i, 0.)) {
        if (phase_i == ChemicalPhase::Fluid) {
          has_fluid = 1;
        } else if (phase_i == ChemicalPhase::Solid) {
          has_solids = 1;
        }
      }
    }

    if (has_fluid && has_solids) {
      m_has_fluid_and_solids = 1;
    } else {
      m_has_fluid_and_solids = 0;
    }
  } else {
    amrex::Abort("Error: unrecognized reaction type");
  }

  return;
}
