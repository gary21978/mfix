#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>

#include <mfix_prop_specific_heat.H>
#include <mfix_reporter.H>
#include <mfix_fix_inputs.H>

using namespace amrex;

int SpecificHeat::
define ( std::string base )
{
  ParmParse ppMFIX("mfix");
  ParmParse ppSpecies("species");
  ParmParse ppBase(base); // e.g., fluid0

  // Check if solving species and if so, get the list defining base
  int solve_species(0);
  ppMFIX.query("solve_species", solve_species);

  Vector<std::string> species;
  if (solve_species) {
    ppBase.queryarr("species", species);
  }

  //
  // base.specific_heat = constant, NASA7-poly, or mixture,
  //
  // |-- base.specific_heat = constant
  // |-- base.specific_heat = NASA-7
  // |-- base.specific_heat = mixture
  //     |
  //     |-- species.specific_heat = constant
  //     +-- species.specific_heat = NASA7-poly
  //
  std::string model_str;
  if (!ppBase.query("specific_heat", model_str)) {
    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Specific heat model not specified for " << base << "\n"
      << "Options: constant, NASA7-poly, mixture\n"
      << "Please correct the input deck.";
    return 1;
  }

  int const is_mixture( toLower(model_str) == "mixture" );
  m_nspecies = (is_mixture) ? static_cast<int>(species.size()) : 0;

  if (is_mixture) {

    if (!solve_species) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Specific heat mixture model requires solving species equations.\n"
        << "Please correct the input deck.";
      return 1;
    }

    if (m_nspecies == 0) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Specific heat mixture model requires solving species equations.\n"
        << "No species found for " << base << "\n"
        << "Please correct the input deck.";
      return 1;
    }

    model_str.clear();
    if (!ppSpecies.query("specific_heat", model_str)) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Specific heat model not specified for species\n"
        << "Options: constant, NASA7-poly\n"
        << "Please correct the input deck.";
      return 1;
    }
  }

  // First read the reference temperature inputs
  m_h_T_offset.resize((amrex::max(m_nspecies,1)) + 1, 0);
  m_d_T_offset.resize((amrex::max(m_nspecies,1)) + 1, 0);

  m_h_coeff_offset.resize((amrex::max(m_nspecies,1)) + 1, 0);
  m_d_coeff_offset.resize((amrex::max(m_nspecies,1)) + 1, 0);

  if (ReadT( base, species ))
  { return 1; }

  // Each NASA-7 polynomial contains 7 coefficients
  // for each temperatures range
  // Note: constant is converted to NASA-7 for use
  m_ncomp = m_h_coeff_offset[(amrex::max(m_nspecies,1))];
  m_d_T.resize(m_h_T.size(), 0.);

  m_h_coeffs.resize(m_ncomp, 0.);
  m_d_coeffs.resize(m_ncomp, 0.);

  if (toLower(model_str) == "constant") {
    m_model = SpecificHeatModel::Constant;

    Real Tref(298.0);
    if (!ppBase.query("reference_temperature", Tref)) {
      reporter::Log(reporter::Info,__FILE__, __LINE__)
        << "Input not provided: " << base << ".reference_temperature\n"
        << "Using default value: " << Tref;
    }

    if (ReadConstant( Tref, base, species ))
    { return 1; }

  } else if (toLower(model_str) == "nasa7-poly") {
    m_model = SpecificHeatModel::NASA7Polynomials;

    if (ReadNASA7( base, species ))
    { return 1; }

  } else {

    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Unknown or invalid specific heat model:" << model_str << "\n"
      << "Options: constant, NASA7-poly\n"
      << "Please correct the input deck.";
    return 1;

  }

  // All the data for specific heat and enthalpy has been read.
  // Perform some checks on the values

  EnthalpyData enthalpy( m_nspecies, m_h_coeffs.data(), m_h_T_offset.data(),
      m_h_coeff_offset.data(), m_h_T.data());
  SpecificHeatData specificHeat( m_nspecies, m_h_coeffs.data(), m_h_T_offset.data(),
      m_h_coeff_offset.data(), m_h_T.data());

  // Sanity checks on specific heat and enthalpy.
  int ignore(0);
  if ( is_mixture ) {
    {
      FixInputs fix("Dec. 2025");
      fix.swap<int>("species.ignore_discontinuities",
                    "species.specific_heat.ignore_discontinuities");
    }
    ppSpecies.query("specific_heat.ignore_discontinuities", ignore);
    for (int n(0); n<m_nspecies; ++n) {
      if (CheckEntry(ignore, n, species[n], enthalpy, specificHeat))
      { return 1; }
    }

  } else {
    {
      FixInputs fix("Dec. 2025");
      fix.swap<int>(base + ".ignore_discontinuities",
                  base + ".specific_heat.ignore_discontinuities");
    }

    ppBase.query("specific_heat.ignore_discontinuities", ignore);
    if (CheckEntry(ignore, 0, base, enthalpy, specificHeat))
    { return 1; }
  }

  // All checks passed, copy from host to device.
  Gpu::copyAsync(Gpu::hostToDevice, m_h_coeffs.begin(), m_h_coeffs.end(), m_d_coeffs.begin());
  Gpu::copyAsync(Gpu::hostToDevice, m_h_T_offset.begin(), m_h_T_offset.end(), m_d_T_offset.begin());
  Gpu::copyAsync(Gpu::hostToDevice, m_h_coeff_offset.begin(), m_h_coeff_offset.end(), m_d_coeff_offset.begin());
  Gpu::copyAsync(Gpu::hostToDevice, m_h_T.begin(), m_h_T.end(), m_d_T.begin());

  return 0;
}


int SpecificHeat::
ReadConstant ( Real const a_Tref, std::string base,
               Vector<std::string> a_species )
{
  if (m_nspecies > 0) {

    for( int n(0); n<m_nspecies; ++n) {

      // species.[species name].specific_heat.constant = <Real>
      std::string species_base( "species."+a_species[n]);
      if (ReadConstantEntry(a_Tref, species_base, n) ) { return 1; }
    }

  } else {

    // [base].specific_heat.constant = <Real>
    if (ReadConstantEntry(a_Tref, base, 0) ) { return 1; }

  }
  return 0;
}


int SpecificHeat::
ReadConstantEntry ( Real const a_Tref, std::string base,
                    int const comp )
{
  ParmParse pp(base);

  // Read constant specific heat
  Real value(0.);
  if (!pp.query("specific_heat.constant", value)) {
    reporter::Log(reporter::Error,__FILE__, __LINE__)
      << "Invalid or missing input: "
      << "Please correct the input deck.";
    return 1;
  }

  int n_T = m_h_T_offset[comp+1] - m_h_T_offset[comp];
  for (int t(0); t <= n_T; ++t) {
    int const cp_coeff(m_h_coeff_offset[comp] + 7*t);

    AMREX_ASSERT(cp_coeff < static_cast<int>(m_h_coeffs.size()));

    m_h_coeffs[cp_coeff] = value;
  }

  // Read enthalpy of formation
  value = 0.;
  if (!pp.query("enthalpy_of_formation", value)) {
    reporter::Log(reporter::Info,__FILE__, __LINE__)
      << "Input not provided: " << base << ".enthalpy_of_formation\n"
      << "Using default value: " << value;
  }

  for (int t(0); t <= n_T; ++t) {
    int const cp_coeff(m_h_coeff_offset[comp] + 7*t + 0);
    int const h_coeff(m_h_coeff_offset[comp] + 7*t + 5);

    AMREX_ASSERT(h_coeff < static_cast<int>(m_h_coeffs.size()));

    m_h_coeffs[h_coeff] = value - a_Tref*m_h_coeffs[cp_coeff];
  }

  return 0;
}

int SpecificHeat::
ReadNASA7 ( std::string base, Vector<std::string> species )
{
  if (m_nspecies) {
    for (int n(0); n<m_nspecies; ++n) {
      std::string species_base( "species."+ species[n]);
      if (ReadNASA7Entry(species_base, n) ) { return 1; }
    }
  } else {
    if (ReadNASA7Entry(base, 0) ) { return 1; }
  }
  return 0;
}


int SpecificHeat::
ReadNASA7Entry ( std::string name, int const comp)
{
  ParmParse pp(name);

  Real MW(0.);
  if (!pp.query("molecular_weight", MW)) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid or missing input: " << name + ".molecular_weight" << "\n"
        << "Please correct the input deck.";
  }

  // specific gas constant (J/kg.K)
  Real const R_MW (8.31446261815324/MW);

  // NASA-7 polynomial:
  // ------------------------------------------------------------------------//
  // Specific Heat:
  // Cp = (a0*T^0 * a1*T^1 + a2*T^2 + a3*T^3 + a4*T^4)*(R/MW)
  //
  // Enthalpy:
  // Hf = (a0*T * a1*T^2/2 + a2*T^3/3 + a3*T^4/4 + a4*T^5/5 + a5)*(R/MW)

  for( int n(0); n<7; ++n) {

    Vector<Real> values;
    std::string NASA7a("specific_heat.NASA7.a" + std::to_string(n));

    if (!pp.queryarr(NASA7a, values)) {

      // The 7-th coefficient is currently not used; therefore, do not
      // throw an error if it was not provided.
      if ( n == 6 ) { continue; }

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid or missing input: " << NASA7a << "\n"
        << "Please correct the input deck.";

      return 1;
    }

    if ( values.size() != (m_h_T_offset[comp+1] - m_h_T_offset[comp] + 1) ) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid or missing input: " << NASA7a << "\n"
        << "Please correct the input deck.";
      return 1;
    }

    for( int t(0); t < values.size(); ++t) {

      int const coeff(m_h_coeff_offset[comp] + 7*t + n);

      AMREX_ASSERT(coeff < static_cast<int>(m_h_coeffs.size()));

      m_h_coeffs[coeff] = values[t]*R_MW;
    }
  }

  return 0;
}

int SpecificHeat::
ReadT ( std::string base, Vector<std::string> a_species )
{
  if (m_nspecies) {

    for( int n(0); n<m_nspecies; ++n) {
      std::string species_base( "species."+a_species[n]);
      if (ReadTEntry(species_base, n) ) { return 1; }
    }

  } else {

    if (ReadTEntry(base, 0) ) { return 1; }

  }
  return 0;
}


int SpecificHeat::
ReadTEntry ( std::string name, int const comp)
{
  ParmParse pp(name);

  Vector<Real> temps;
  if (!pp.queryarr("specific_heat.NASA7.Tsplit", temps)) {
    m_h_T_offset[comp + 1] = m_h_T_offset[comp] + 1;
    m_h_coeff_offset[comp + 1] = m_h_coeff_offset[comp] + 14;
    m_h_T.push_back(1000.);  // default temp of 1000
  } else if ((temps.size() == 1) && (temps[0] == -1.)) {
    m_h_T_offset[comp + 1] = m_h_T_offset[comp];
    m_h_coeff_offset[comp + 1] = m_h_coeff_offset[comp] + 7;
  } else {
    m_h_T_offset[comp + 1] = m_h_T_offset[comp] + static_cast<int>(temps.size());
    m_h_coeff_offset[comp + 1] = m_h_coeff_offset[comp] + 7*(1 + static_cast<int>(temps.size()));
    for (int t(0); t < temps.size() ; ++t) {
      m_h_T.push_back(temps[t]);
    }
  }
  return 0;
}

int SpecificHeat::
CheckEntry ( int const a_ignore, int const comp, std::string base,
              EnthalpyData& a_H, SpecificHeatData& a_Cp)
{
  int logLevel = ((a_ignore == 1) ? reporter::Warning : reporter::Error);

  int n_T = a_Cp.p_T_offset[comp+1] - a_Cp.p_T_offset[comp];

  for (int t(0); t < n_T; ++t) {
    const Real loT = (1. - std::numeric_limits<amrex::Real>::epsilon()) * a_Cp.p_T[a_Cp.p_T_offset[comp] + t];
    const Real hiT = (1. + std::numeric_limits<amrex::Real>::epsilon()) * a_Cp.p_T[a_Cp.p_T_offset[comp] + t];

    { Real const loCp = a_Cp(comp, loT);
      Real const hiCp = a_Cp(comp, hiT);

      Real const max_abs_Cp = amrex::max(std::abs(loCp), std::abs(hiCp));
      Real const rel_error_Cp = std::abs((hiCp - loCp) / max_abs_Cp);

      if (rel_error_Cp > 1.e-6) {
        reporter::Log(logLevel, __FILE__, __LINE__)
          << "Specific heat is not continuous at T=" << a_Cp.p_T[a_Cp.p_T_offset[comp]] << " for " << base << "\n"
          << "The relative error is " << rel_error_Cp*100 << "%\n"
          << "Please correct the input deck.";
          if (!a_ignore) { return 1;}
      }
    }

    { Real const loH = a_H(comp, loT);
      Real const hiH = a_H(comp, hiT);

      Real const max_abs_H = amrex::max(std::abs(loH), std::abs(hiH));
      Real const rel_error_H = std::abs((hiH - loH) / max_abs_H);

      if (rel_error_H > 1.e-6) {
        reporter::Log(logLevel, __FILE__, __LINE__)
          << "Enthapy heat is not continuous at T=" << a_Cp.p_T[a_Cp.p_T_offset[comp]] << " for " << base << "\n"
          << "The relative error is " << rel_error_H*100 << "%\n"
          << "Please correct the input deck.";
          if (!a_ignore) { return 1;}
      }
    }
  }

  return 0;
}
