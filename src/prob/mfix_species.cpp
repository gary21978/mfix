#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <mfix_species.H>
#include <mfix_fluid.H>
#include <mfix_reporter.H>


using namespace amrex;


MFIXSpecies::MFIXSpecies()
  : m_diffusivity_model(DiffusivityModel::Undefined)
  , m_solve(0)
  , m_nspecies(0)
  , m_names(0)
  , m_IDs(0)
  , m_MW_k(0)
  , m_D(0)
  , m_is_initialized(0)
{}


void
MFIXSpecies::Initialize ()
{
  // Set the initialization flag
  m_is_initialized = 1;

  int solve_enthalpy(0);

  amrex::ParmParse ppMFIX("mfix");
  ppMFIX.query("advect_enthalpy", solve_enthalpy);

  amrex::ParmParse pp("species");

  if (pp.contains("solve")) {

    // Get the species names and store them in an alphabetical ordered way
    pp.getarr("solve", m_names);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_names.size() > 0,
                                     "No input provided for species.solve");

    // Disable the species solver if the species are defined as "None" (case
    // insensitive) or 0
    if (amrex::toLower(m_names[0]).compare("none") == 0) {
      m_solve = 0;

    } else {
      m_solve = 1;
      m_nspecies = m_names.size();

      m_IDs.resize(m_nspecies);
      for (int n(0); n < m_nspecies; n++) {
        m_IDs[n] = n;
      }

      m_MW_k.resize(m_nspecies, 0.);
      m_densities_k.resize(m_nspecies, 0.);

    }

    if (m_solve) {
      bool contains_densities(false);

      // Get molecular weights input --------------------------------//
      for (int n(0); n < m_nspecies; n++) {
        std::string species_prefix = "species." + m_names[n];
        amrex::ParmParse ppSpecies(species_prefix);

        bool contains_MW_k = ppSpecies.query("molecular_weight", m_MW_k[n]);

        if (!contains_MW_k) {
          if (amrex::ParallelDescriptor::IOProcessor()) {
            std::string message = "Input not provided. Assuming MW_" + m_names[n] + " = 0";
            amrex::Warning(message);
          }
        }

        bool contains_density_k = ppSpecies.query("density", m_densities_k[n]);
        contains_densities = contains_densities || contains_density_k;
      }

      if (!contains_densities) {
        m_densities_k.clear();
        m_densities_k = {};
      }

      // Get diffusivity model input --------------------------------//
      {
        std::string model;
        pp.get("diffusivity", model);

        if (amrex::toLower(model) == "constant") {
          m_diffusivity_model = DiffusivityModel::Constant;

          pp.get("diffusivity.constant", m_D);

        } else {
          amrex::Abort("Unknown species mass diffusivity model");
        }
      }

    } // solve species
  }
}
