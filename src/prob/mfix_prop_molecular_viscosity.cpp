#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>

#include <mfix_prop_molecular_viscosity.H>
#include <mfix_reporter.H>

using namespace amrex;

int MolViscosity::
define ( const std::string& fluid_name)
{

  ParmParse ppMFIX("mfix");
  ParmParse ppSpecies("species");
  ParmParse ppFluid(fluid_name); // e.g., fluid0

  // Check if solving a mixture viscosity
  std::string model;
  ppFluid.query("viscosity.molecular", model);
  bool is_mixture = (toLower(model) == "mixture");


  Vector<std::string> species;
  if (is_mixture) {
    // Check if solving species and if so, get the list defining base
    int solve_species(0);
    ppMFIX.query("solve_species", solve_species);

    if (solve_species) {
      ppFluid.queryarr("species", species);
    } else {
       reporter::Log(reporter::Error,__FILE__, __LINE__)
         << "Mixture fluid viscosity model requires the fluid to be a mixture.\n"
         << "Please correct the input deck.";
       return 1;
    }
  }

  m_ncomp = amrex::max(1, static_cast<int>(species.size()));

  // Get the specified molecular viscosity model
  {
    ParmParse ppBase = (m_ncomp > 1) ? ppSpecies : ppFluid;

    if (ppBase.query("viscosity.molecular", model)){
      std::string mlower = toLower(model);
      if (mlower == "none") {
        m_model = MolViscosityModel::Undefined;
      } else if (mlower == "constant") {
        m_model = MolViscosityModel::Constant;
        m_h_coeffs.resize(m_ncomp);
        m_d_coeffs.resize(m_ncomp);
      } else if (mlower == "sutherland") {
        m_model = MolViscosityModel::Sutherland;
        m_h_coeffs.resize(m_ncomp*2);
        m_d_coeffs.resize(m_ncomp*2);
      } else if (mlower == "reid") {
        m_model = MolViscosityModel::Reid;
        m_h_coeffs.resize(m_ncomp*4);
        m_d_coeffs.resize(m_ncomp*4);
      } else {
        amrex::Abort("Unknown molecular viscosity model " + model);
      }
    } else {
      m_model = MolViscosityModel::Undefined;
    }
  }


  // Get specified parameters for molecular viscosity
  std::string prefix = (species.size() > 0) ? "species." : "";
  const auto& names = (species.size() > 0) ? species : Vector<std::string>{fluid_name};
  for (int n(0); n < m_ncomp; ++n) {
    std::string name = prefix + names[n];

    ParmParse ppBase(name);

    if (m_model == MolViscosityModel::Constant) {
      if (!ppBase.query("viscosity.molecular.constant", m_h_coeffs[n])) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Constant molecular viscosity is not defined for "
          << names[n] << ". Please correct the input deck.";
        return 1;
      }
    } else if  (m_model == MolViscosityModel::Sutherland) {
      amrex::Real T_ref = 0.;
      amrex::Real mu_ref = 0.;
      m_h_coeffs[n*2 + 1] = 0.; // Sutherland temperature

      if (!ppBase.query("viscosity.molecular.Sutherland.T_ref", T_ref) ) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Sutherland reference temperature is not defined for "
          << names[n] << ". Please correct the input deck.";
        return 1;
      }

      if (!ppBase.query("viscosity.molecular.Sutherland.mu_ref", mu_ref) ) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Sutherland reference viscosity is not defined for "
          << names[n] << ". Please correct the input deck.";
        return 1;
      }

      if (!ppBase.query("viscosity.molecular.Sutherland.S", m_h_coeffs[n*2 + 1]) ) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Sutherland temperature is not defined for "
          << names[n] << ". Please correct the input deck.";
        return 1;
      }

      m_h_coeffs[n*2] = (mu_ref * (T_ref + m_h_coeffs[n*2 + 1])) /
         std::pow(T_ref, (3./2.)); // Sutherland constant

    } else if (m_model == MolViscosityModel::Reid) {
      if (!ppBase.query("viscosity.molecular.Reid.A", m_h_coeffs[n*4]) ) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Reid viscosity constant A is not defined for "
          << names[n] << ". Please correct the input deck.";
        return 1;
      }

      if (!ppBase.query("viscosity.molecular.Reid.B", m_h_coeffs[n*4 + 1]) ) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Reid viscosity constant B is not defined for "
          << names[n] << ". Please correct the input deck.";
        return 1;
      }

      if (!ppBase.query("viscosity.molecular.Reid.C", m_h_coeffs[n*4 + 2]) ) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Reid viscosity constant C is not defined for "
          << names[n] + ". Please correct the input deck.";
        return 1;
      }

      if (!ppBase.query("viscosity.molecular.Reid.D", m_h_coeffs[n*4 + 3]) ) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Reid viscosity constant D is not defined for "
          << names[n] << ". Please correct the input deck.";
        return 1;
      }
    } else if (m_model == MolViscosityModel::Usr) {
      // Read user species molecular viscosity inputs here
    } else {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Viscosity not defined for "
        << names[n] << ". Please correct the input deck.";
      return 1;
    }

  }

  m_defined = 1;
  Gpu::copyAsync(Gpu::hostToDevice, m_h_coeffs.begin(),
      m_h_coeffs.end(), m_d_coeffs.begin());

  return 0;
}
