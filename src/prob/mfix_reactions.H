#ifndef MFIX_REACTIONS_H_
#define MFIX_REACTIONS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Gpu.H>

#include <mfix_species.H>


struct ReactionType { enum {Undefined=-1, Eulerian, Lagrangian}; };


struct ChemicalPhase { enum {Undefined=-1, Fluid, Solid}; };


struct InterphaseChemTxfrIndexes {

  InterphaseChemTxfrIndexes (const int nspecies_g,
                             const int solve_reactions)
    : chem_ro_gk(0)
    , chem_h(chem_ro_gk + nspecies_g*(solve_reactions>0))
    , count(chem_h + 1*(solve_reactions>0)) // Num of chem rxns txfr terms
  {}

  const int chem_ro_gk;
  const int chem_h;
  const int count;
};


// Forward declaration
class MFIXChemicalReaction;


namespace reactions_algorithms {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int linear_search (const char* names_array,
                   const int* names_idxs,
                   const int size,
                   const char* name)
{
  int name_size(0);
  while (name[name_size] != '\0')
  {
    AMREX_ASSERT(name_size < 80);
    name_size++;
  }
  AMREX_ASSERT(name_size > 0);

  for (int i(0); i < size; ++i) {

    const char* name_string = &(names_array[names_idxs[i]]);
    const int length = names_idxs[i+1] - names_idxs[i];

    if (length == name_size) {
      int match = 1;
      for (int n(0); n < length; ++n) {
        if (name_string[n] != name[n]) {
          match = 0;
          break;
        }
      }

      if (match)
        return i;
    }
  }

  return -1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int binary_search (const char* names_array,
                   const int* names_idxs,
                   const int size,
                   const char* name)
{
  int left = 0;
  int right = size-1;

  int name_size(0);
  while (name[name_size] != '\0')
  {
    AMREX_ASSERT(name_size < 80);
    name_size++;
  }
  AMREX_ASSERT(name_size > 0);

  while (left <= right) {
    const int mid = left + (right - left) / 2;

    const char* name_string = &(names_array[names_idxs[mid]]);
    const int length = names_idxs[mid+1] - names_idxs[mid];

    if (length == name_size) {
      int match = 1;
      for (int n(0); n < length; ++n) {
        if (name_string[n] != name[n]) {
          match = 0;
          break;
        }
      }

      if (match)
        return mid;
    }

    {
      for (int n(0); n < amrex::max(length, name_size); ++n) {
        if (n >= length) {
          left = mid+1;
          break;
        } else if (n >= name_size) {
          right = mid-1;
          break;
        }

        if (name_string[n] == name[n]) {
          continue;
        } else if (name_string[n] < name[n]) {
          left = mid+1;
          break;
        } else {
          right = mid-1;
          break;
        }
      }
    }
  }

  return -1;
}

} // end namespace reactions_algorithms


class MFIXReactionsParms {

  public:
    AMREX_GPU_HOST_DEVICE
    MFIXReactionsParms ()
      : m_nreactions(0)
    {}

    AMREX_GPU_HOST_DEVICE
    MFIXReactionsParms (const int nreactions)
      : m_nreactions(nreactions)
    {}

    int m_nreactions;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int reaction_idx (const char* reaction_name) const
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXReactionsParms::reaction_idx()");
#endif

      int ordered_idx = reactions_algorithms::linear_search(m_reactions_names,
          m_reactions_names_idxs, m_nreactions, reaction_name);
      return m_reactions_names_IDs[ordered_idx];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int has_fluid_and_solids (const int q) const
    {
      AMREX_ASSERT(q >= 0 && q < m_nreactions);
      return m_reactions_have_fluid_and_solids[q];
    }

  private:
    // Friendship so MFIXFluidPhase can access private members
    friend class MFIXReactions;

    char* m_reactions_names;
    int*  m_reactions_names_idxs;
    int*  m_reactions_names_IDs;
    int*  m_reactions_have_fluid_and_solids;
};


class MFIXReactions {

  public:
    // Constructor
    MFIXReactions ();

    // Destructor
    ~MFIXReactions ();

    // Initialization procedure
    void Initialize (const MFIXSpecies& species);

    int solve () const { return m_solve; }

    int lagrangian_nreactions () const { return m_lagrangian_nreactions; }

    int eulerian_nreactions () const { return m_eulerian_nreactions; }

    const std::vector<std::string>& reactions () const { return m_reactions_names; }

    const std::vector<std::string>& reaction_equations () const { return m_reactions_equations; }

    int isInitialized () const { return m_is_initialized; }

    template <amrex::RunOn run_on>
    const MFIXReactionsParms& eulerian_parameters () const
    {
      if constexpr (run_on == amrex::RunOn::Cpu)
        return *m_h_eulerian_parameters;
      else if constexpr (run_on == amrex::RunOn::Gpu)
        return *m_d_eulerian_parameters;
    }

    template <amrex::RunOn run_on>
    const MFIXReactionsParms& lagrangian_parameters () const
    {
      if constexpr (run_on == amrex::RunOn::Cpu)
        return *m_h_lagrangian_parameters;
      else if constexpr (run_on == amrex::RunOn::Gpu)
        return *m_d_lagrangian_parameters;
    }

    // Getter
    MFIXChemicalReaction* get_lagrangian (const int q) const
    {
      AMREX_ALWAYS_ASSERT(q < m_lagrangian_nreactions);
      return m_lagrangian_chemical_reactions[q];
    }

    // Getter
    MFIXChemicalReaction* get_eulerian (const int q) const
    {
      AMREX_ALWAYS_ASSERT(q < m_eulerian_nreactions);
      return m_eulerian_chemical_reactions[q];
    }

  private:
    // Switch for turning on/off chemical reactions modeling
    int m_solve;

    // Number of chemical reactions allowed by the model
    int m_lagrangian_nreactions;
    int m_eulerian_nreactions;

    // Names of chemical reactions allowed by the model
    amrex::Vector<std::string> m_reactions_names;

    // Chemical reactions equations
    amrex::Vector<std::string> m_reactions_equations;

    amrex::Gpu::HostVector<char> m_h_eulerian_reactions_names;
    amrex::Gpu::HostVector<int> m_h_eulerian_reactions_names_idxs;
    amrex::Gpu::HostVector<int> m_h_eulerian_reactions_names_IDs;
    amrex::Gpu::DeviceVector<char> m_d_eulerian_reactions_names;
    amrex::Gpu::DeviceVector<int> m_d_eulerian_reactions_names_idxs;
    amrex::Gpu::DeviceVector<int> m_d_eulerian_reactions_names_IDs;

    amrex::Gpu::HostVector<char> m_h_lagrangian_reactions_names;
    amrex::Gpu::HostVector<int> m_h_lagrangian_reactions_names_idxs;
    amrex::Gpu::HostVector<int> m_h_lagrangian_reactions_names_IDs;
    amrex::Gpu::HostVector<int> m_h_lagrangian_reactions_have_fluid_and_solids;
    amrex::Gpu::DeviceVector<char> m_d_lagrangian_reactions_names;
    amrex::Gpu::DeviceVector<int> m_d_lagrangian_reactions_names_idxs;
    amrex::Gpu::DeviceVector<int> m_d_lagrangian_reactions_names_IDs;
    amrex::Gpu::DeviceVector<int> m_d_lagrangian_reactions_have_fluid_and_solids;

    // Reactions unique identifying code
    amrex::Vector<int> m_reactions_IDs;

    // Total number of reactions
    int m_nreactions;

    int m_is_initialized;

    MFIXReactionsParms* m_h_eulerian_parameters;
    MFIXReactionsParms* m_d_eulerian_parameters;

    MFIXReactionsParms* m_h_lagrangian_parameters;
    MFIXReactionsParms* m_d_lagrangian_parameters;

    // Vector of eulerian chemical reactions
    amrex::Vector<MFIXChemicalReaction*> m_eulerian_chemical_reactions;

    // Vector of lagrangian chemical reactions
    amrex::Vector<MFIXChemicalReaction*> m_lagrangian_chemical_reactions;
};


// Class for chemical reaction details
class MFIXChemicalReaction {

  public:
    // Default constructor
    explicit MFIXChemicalReaction () = default;

    // class constructor
    explicit MFIXChemicalReaction (const std::string& name,
                                   const std::string& reaction,
                                   const MFIXSpecies& species);

    MFIXChemicalReaction (const MFIXChemicalReaction&) = delete;

    MFIXChemicalReaction& operator = (const MFIXChemicalReaction&) = delete;

    // Default destructor
    ~MFIXChemicalReaction () = default;

    int get_type () { return m_type; }

    std::string const& get_name () const { return m_name; }

    const amrex::Vector<int>& get_phases () const { return m_phases; }

    const amrex::Vector<std::string>& get_reactants () const { return m_reactants; }

    const amrex::Vector<int>& get_reactants_ids () const { return m_reactants_IDs; }

    const amrex::Vector<amrex::Real>& get_reactants_coeffs () const { return m_reactants_coeffs; }

    const amrex::Vector<int>& get_reactants_phases () const { return m_reactants_phases; }

    const amrex::Vector<std::string>& get_products () const { return m_products; }

    const amrex::Vector<int>& get_products_ids () const { return m_products_IDs; }

    const amrex::Vector<amrex::Real>& get_products_coeffs () const { return m_products_coeffs; }

    const amrex::Vector<int>& get_products_phases () const { return m_products_phases; }

    int has_fluid_and_solids () const { return m_has_fluid_and_solids; }

  private:
    // Reaction type [eulerian/lagrangian]
    int m_type;
    int m_has_fluid_and_solids;

    // chemical formula
    const std::string m_name;
    const std::string m_formula;

    // Phases
    amrex::Vector<int> m_phases;

    // Reactants
    amrex::Vector<std::string> m_reactants;
    amrex::Vector<int> m_reactants_IDs;
    amrex::Vector<amrex::Real> m_reactants_coeffs;
    amrex::Vector<int> m_reactants_phases;

    // Products
    amrex::Vector<std::string> m_products;
    amrex::Vector<int> m_products_IDs;
    amrex::Vector<amrex::Real> m_products_coeffs;
    amrex::Vector<int> m_products_phases;

    // Get reactants data
    std::string parse_reactants (const std::string& formula);

    // Get reactants data
    std::string parse_products (const std::string& formula);

    // Extract stoichiometric data
    void parse_stoichiometric_data (const std::string& s,
                                    amrex::Vector<std::string>& compounds,
                                    amrex::Vector<int>& compounds_id,
                                    amrex::Vector<amrex::Real>& coefficients,
                                    amrex::Vector<int>& phases,
                                    const MFIXSpecies& species);

    void check_mass_balance (const MFIXSpecies& species);

    // Parse reaction
    void parse_reaction (const MFIXSpecies& species);

    amrex::Real m_mass_balance_tolerance;
};

#endif
