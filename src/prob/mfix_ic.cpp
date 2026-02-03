#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>

#include <AMReX_ParmParse.H>

#include <mfix_reporter.H>
#include <mfix_ic.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_regions.H>
#include <mfix_species.H>

#include <algorithm>


using namespace amrex;


void
MFIXInitialConditions::
Initialize (const MFIXRegions& regions,
            MFIXFluidPhase& fluid,
            MFIXSolidsPhase& solids,
            MFIXDEM& dem,
            MFIXPIC& pic)
{
  // The default type is "ASCIIFile" but we can overwrite that in the inputs
  // file with "Auto"

  amrex::ParmParse ppMFIX("mfix");

  std::string init_type = "ASCIIFile";
  ppMFIX.query("particle_init_type", init_type);

  if (toLower(init_type) == "asciifile") {
    m_particle_init_type = ParticleInitType::AsciiFile;

  } else if (toLower(init_type) == "auto") {
    m_particle_init_type = ParticleInitType::Auto;

  } else {
     amrex::Abort("Bad particle_init_type " + init_type);
  }

  amrex::ParmParse pp("ic");

  std::vector<std::string> input_regions;
  pp.queryarr("regions", input_regions);

  {
    amrex::ParmParse ppIC("ic");

    ppIC.query("allow_regions_overlap", m_allow_overlap);

    std::string ranking_type("Inputs");
    ppIC.query("ranking_type", ranking_type) ;
    std::string rlower(toLower(ranking_type));
    if(rlower == "inputs") {
      m_ic_ranking.set_type(ICRankingType::Inputs);
    } else if (rlower == "volume") {
      m_ic_ranking.set_type(ICRankingType::Volume);
    } else if (rlower == "priority") {
      m_ic_ranking.set_type(ICRankingType::Priority);
    } else {
      amrex::Abort("Unknown IC ranking type " + ranking_type);
    }
  }

  // Loop over ICs
  for (size_t icv=0; icv < input_regions.size(); icv++) {

    Real volfrac_total(0.0);

    IC_t new_ic(&fluid);

    // Set the region for the initial condition.
    new_ic.region = regions.getRegion(input_regions[icv]);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(new_ic.region != nullptr, "Invalid IC region");

    // Check IC regions ranking inputs
    {
      std::string ic_region_prefix = "ic."+input_regions[icv];
      amrex::ParmParse ppICRegion(ic_region_prefix);

      if(m_ic_ranking.get_type() == ICRankingType::Inputs) {
        new_ic.inputs_order = icv;
      } else if (m_ic_ranking.get_type() == ICRankingType::Volume) {
        new_ic.volume = new_ic.region->volume();
      } else if (m_ic_ranking.get_type() == ICRankingType::Priority) {
        ppICRegion.query("priority", new_ic.priority);
      }

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(new_ic.priority >= 1, "IC priority must be >= 1");
    }

    // Get fluid data.
    if (fluid.solve()) {

      std::string ic_region_fluid_prefix = "ic."+input_regions[icv]+"."+fluid.names(0);
      amrex::ParmParse ppFluid(ic_region_fluid_prefix);

      ppFluid.get("volfrac", new_ic.fluid.volfrac);
      volfrac_total += new_ic.fluid.volfrac;

      if (new_ic.fluid.velocity.get(ppFluid, "velocity")) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Failed to process initial condition fluid velocity.\n"
          << "IC region: " << input_regions[icv] << '\n'
          << "Please correct the input deck.";
      }

      if (fluid.solve_species()) {

        const int nspecies_g = fluid.nspecies();
        new_ic.fluid.species.resize(nspecies_g, FVAR_({1}));

        std::string ic_region_fluid_species_prefix = ic_region_fluid_prefix+".species";
        amrex::ParmParse ppSpecies(ic_region_fluid_species_prefix);

        // Auxiliary variable to check that species sum up to 1
        Real total_mass_fraction(0);

        for (int n(0); n < nspecies_g; n++) {

          std::string fluid_species = fluid.species_names(n);
          if (new_ic.fluid.species[n].get(ppSpecies, fluid_species)) {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Failed to process initial condition fluid species.\n"
              << "IC region: " << input_regions[icv] << "   "
              << "species: " << fluid_species << '\n'
              << "Please correct the input deck.";
          }

          total_mass_fraction += new_ic.fluid.get_species(n);
        }

        // Sanity check that the input species mass fractions sum up to 1
        if ( !amrex::almostEqual(total_mass_fraction,1.) ) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Fluid species mass fractions do not sum to 1.\n"
            << "IC region: " << input_regions[icv] << "   "
            << "sum(X): " << total_mass_fraction << '\n'
            << "Please correct the input deck.";
        }
      }

      // Get density
      if( new_ic.fluid.density.get(ppFluid, "density")) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Failed to process initial condition fluid density.\n"
          << "IC region: " << input_regions[icv] << '\n'
          << "Please correct the input deck.";
      }

      // Get temperature
      if( new_ic.fluid.temperature.get(ppFluid, "temperature")) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Failed to process initial condition fluid temperature.\n"
          << "IC region: " << input_regions[icv] << '\n'
          << "Please correct the input deck.";
      }

      // Get pressure
      new_ic.fluid.pressure_defined = ppFluid.query("pressure", new_ic.fluid.pressure);

      if (fluid.constraint.isIncompressibleFluid() ) {

        if (!new_ic.fluid.density.is_defined()) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "ICs must include density with IncompressibleFluid constraint\n"
            << "IC region: " << input_regions[icv] << '\n'
            << "Please correct the input deck.";
        }

      } else {

        if (new_ic.fluid.density.is_defined()) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "ICs cannot define density with IdealGas constraints\n"
            << "IC region: " << input_regions[icv] << '\n'
            << "Please correct the input deck.";
        }

        if (!fluid.p_therm_defined()) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Thermodynamic pressure is required for IdealGas constraints\n"
            << "Please correct the input deck.";
        }
      }
    }

    amrex::Real gran_temp(0.);

    if (dem.solve() || pic.solve()) {

      // solids volfrac in IC region is > 0
      if (new_ic.fluid.volfrac < 1.0) {

        // Get the list of solids used in defining the IC region
        std::vector<std::string> solids_names;
        {
          std::string ic_region_prefix = "ic."+input_regions[icv];
          amrex::ParmParse ppICRegion(ic_region_prefix);
          ppICRegion.getarr("solids", solids_names);

          if (solids_names.size() != 1) {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Initial conditions only support one solid per region.\n"
              << "IC region: " << input_regions[icv] << "   "
              << "solids: " << solids_names.size() << '\n'
              << "Please correct the input deck.";
          }

          if(AutoParticleInit()) {
            ppICRegion.get("packing", new_ic.packing);
            ppICRegion.query("granular_temperature", gran_temp);
          }
        }

        for (size_t lcs(0); lcs < solids_names.size(); ++lcs) {

          SOLIDS_t new_solid;

          std::string ic_region_solid_prefix = "ic."+input_regions[icv]+"."+solids_names[lcs];
          amrex::ParmParse PPICRegionSolid(ic_region_solid_prefix);

          new_solid.name = solids_names[lcs];
          new_solid.phase = solids.name_to_phase(solids_names[lcs]);

          if (solids.solve_species())
          {
            std::string ic_region_solid_species_prefix = ic_region_solid_prefix + ".species";
            amrex::ParmParse ppSpecies(ic_region_solid_species_prefix);

            new_solid.species.resize(solids.nspecies());

            amrex::Real total_mass_fraction(0);

            for (int n(0); n < solids.nspecies(); n++) {

              new_solid.species[n].mass_fraction = 0.;

              std::string current_species = solids.species_names(n);
              ppSpecies.query(current_species, new_solid.species[n].mass_fraction);

              total_mass_fraction += new_solid.species[n].mass_fraction;
            }

            // Sanity check that the input species mass fractions sum up to 1
            if ( !amrex::almostEqual(total_mass_fraction,1.0) ) {
                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "Solids mass fractions do not sum to 1.\n"
                  << "IC region: " << input_regions[icv] << "   "
                  << "solids: " << solids_names[lcs] << '\n'
                  << "Please correct the input deck.";
            }

          }

          //
          if(AutoParticleInit()) {

            { // Get solids volume fraction.
              if(!PPICRegionSolid.query("volfrac", new_solid.volfrac)) {
                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "Volume fraction not specified.\n"
                  << "IC region: " << input_regions[icv] << "   "
                  << "solids: " << solids_names[lcs] << '\n'
                  << "Please correct the input deck.";
              }
              volfrac_total += new_solid.volfrac;
            }

            { // Get solids velocity.
              if(!PPICRegionSolid.queryarr("velocity", new_solid.velocity, 0, 3)) {
                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "Solids velocity not specified.\n"
                  << "IC region: " << input_regions[icv] << "   "
                  << "solids: " << solids_names[lcs] << '\n'
                  << "Please correct the input deck.";
              }
            }

            { // Get information about diameter distribution.
              if(new_solid.diameter.define(ic_region_solid_prefix, "diameter",
                   new_solid.get_diameter_data(), dem.solve())) {
                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "Invalid particle diameter parameters.\n"
                  << "IC region: " << input_regions[icv] << "   "
                  << "solids: " << solids_names[lcs] << '\n'
                  << "Please correct the input deck.";
              }
            }

            { // Get information about density distribution
              if (solids.density_model() == MFIXSolidsPhase::DensityModel::SpeciesMixture)
              {
                amrex::Real p_density_inv(0.);
                for (int n(0); n < solids.nspecies(); n++) {
                  p_density_inv += new_solid.species[n].mass_fraction / solids.species_density(n);
                }
                amrex::Real p_density = 1./p_density_inv;

                std::string density_type("");
                const bool contains_density = PPICRegionSolid.query("density", density_type);

                if (contains_density) {
                  if (amrex::toLower(density_type) == "constant") {
                    amrex::Real input_density(0.);
                    PPICRegionSolid.get("density.constant", input_density);

                    if (!amrex::almostEqual(p_density, input_density)) {
                      reporter::Log(reporter::Error,__FILE__, __LINE__)
                        << "Invalid particle density parameters.\n"
                        << "Input particle density "
                        << std::scientific
                        << std::setprecision(std::numeric_limits<double>::max_digits10 - 1)
                        << input_density
                        << "\ndiffers from computed mixture density "
                        << std::scientific
                        << std::setprecision(std::numeric_limits<double>::max_digits10 - 1)
                        << p_density << "\n"
                        << "IC region: " << input_regions[icv] << "   "
                        << "solids: " << solids_names[lcs] << '\n'
                        << "Please correct the input deck.";
                    }
                  } else {
                    reporter::Log(reporter::Error,__FILE__, __LINE__)
                      << "Invalid particle density parameters.\n"
                      << "When solids density model is a mixture of species densities\n"
                      << "the IC solids density should either not be specified"
                      << " or specified as constant.\n"
                      << "IC region: " << input_regions[icv] << "   "
                      << "solids: " << solids_names[lcs] << '\n'
                      << "Please correct the input deck.";
                  }
                } else {
                  PPICRegionSolid.add("density", std::string("constant"));
                  PPICRegionSolid.add("density.constant", amrex::Real(p_density));
                }
              }

              if( new_solid.density.define(ic_region_solid_prefix, "density") ) {
                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "Invalid particle density parameters.\n"
                  << "IC region: " << input_regions[icv] << "   "
                  << "solids: " << solids_names[lcs] << '\n'
                  << "Please correct the input deck.";
              }
            }
          }

          if (fluid.solve_enthalpy() || solids.solve_enthalpy()) {
            PPICRegionSolid.query("temperature", new_solid.temperature);
          }

          new_ic.solids.push_back(new_solid);

        }
      }
      // either fluid.volfrac == 1, or we read particles from particle_input.dat
      else {

        // Get the list of solids used in defining the IC region
        std::vector<std::string> solids_names(0);
        {
          std::string ic_region_prefix = "ic."+input_regions[icv];
          amrex::ParmParse ppICRegion(ic_region_prefix);
          ppICRegion.queryarr("solids", solids_names);

          if (solids_names.size() != 0 && solids_names.size() != 1) {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Initial conditions only support one solid per region.\n"
              << "IC region: " << input_regions[icv] << "   "
              << "solids: " << solids_names.size() << '\n'
              << "Please correct the input deck.";
          }
        }

        for(size_t lcs(0); lcs < solids_names.size(); ++ lcs) {

          SOLIDS_t new_solid;

          std::string ic_region_solid_prefix = "ic."+input_regions[icv]+"."+solids_names[lcs];
          amrex::ParmParse ppICRegionSolid(ic_region_solid_prefix);

          new_solid.name = solids_names[lcs];
          new_solid.phase = solids.name_to_phase(solids_names[lcs]);

          if (fluid.solve_enthalpy() || solids.solve_enthalpy()) {
            ppICRegionSolid.get("temperature", new_solid.temperature);
          }

          if (solids.solve_species()) {

            std::string ic_region_solid_species_prefix = ic_region_solid_prefix+".species";
            amrex::ParmParse ppICRegionSolidSpecies(ic_region_solid_species_prefix);

            new_solid.species.resize(solids.nspecies());

            for (int n(0); n < solids.nspecies(); n++) {
              std::string current_species = solids.species_names(n);
              ppICRegionSolidSpecies.query(current_species,
                                           new_solid.species[n].mass_fraction);
            }
          }

          new_ic.solids.push_back(new_solid);
        }
      }
    }

    if (fluid.solve() && (dem.solve() || pic.solve())) {
      if (AutoParticleInit() && !almostEqual(volfrac_total,1.0)) {
            reporter::Log(reporter::Warning)
              << "Fluid and particle volume fractions do not sum to one.\n"
              << "IC region: " << input_regions[icv]
              << "  Total volume fraction: " << volfrac_total;
      }
    }

    m_ic.push_back(new_ic);
    m_np.push_back(0);
    m_granular_temperature.push_back(gran_temp);

  }

  // Sort the IC on the basis of inputs order, priority or volume
  std::sort(m_ic.begin(), m_ic.end(), m_ic_ranking);

  // This completes the copyAsync to Gpu::DeviceVectors and
  // sets pointers for runtime use.
  for (int icv(0); icv<m_ic.size(); ++icv) {
    for (int m(0); m<m_ic[icv].solids.size(); ++m) {
      m_ic[icv].solids[m].copyAsync();
    }
  }
}


bool
MFIXInitialConditions::ICRanking::operator() (const IC_t& left,
                                              const IC_t& right) const
{
  AMREX_ALWAYS_ASSERT(m_ranking_type != ICRankingType::Undefined);

  if (m_ranking_type == ICRankingType::Inputs) {
    return left.inputs_order < right.inputs_order;
  }

  if (m_ranking_type == ICRankingType::Volume) {
    return left.volume < right.volume;
  }

  if (m_ranking_type == ICRankingType::Priority) {
    return left.priority < right.priority;
  }

  return false;
}
