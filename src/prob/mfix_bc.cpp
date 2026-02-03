#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_Geometry.H>

#include <AMReX_ParmParse.H>

#include <mfix_reporter.H>
#include <mfix_bc.H>
#include <mfix_eb_parms.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_regions.H>
#include <mfix_species.H>

using namespace amrex;


// Constructor
BCList::BCList (const int a_max_levels)
{
  bc_ilo.resize(a_max_levels, nullptr);
  bc_ihi.resize(a_max_levels, nullptr);
  bc_jlo.resize(a_max_levels, nullptr);
  bc_jhi.resize(a_max_levels, nullptr);
  bc_klo.resize(a_max_levels, nullptr);
  bc_khi.resize(a_max_levels, nullptr);
}


// Destructor
BCList::~BCList()
{
  for (auto& lev_bc : bc_ilo) { delete lev_bc; }
  for (auto& lev_bc : bc_ihi) { delete lev_bc; }
  for (auto& lev_bc : bc_jlo) { delete lev_bc; }
  for (auto& lev_bc : bc_jhi) { delete lev_bc; }
  for (auto& lev_bc : bc_klo) { delete lev_bc; }
  for (auto& lev_bc : bc_khi) { delete lev_bc; }
}


void BCList::
MakeBCArrays ( int const a_lev, int const nghost,
               const Geometry& a_geom, bool const a_debug)
{
  if (a_debug) { Print() << "BCList::MakeBCArrays\n"; }

  if (bc_ilo[a_lev] != nullptr) { delete bc_ilo[a_lev]; }
  if (bc_ihi[a_lev] != nullptr) { delete bc_ihi[a_lev]; }
  if (bc_jlo[a_lev] != nullptr) { delete bc_jlo[a_lev]; }
  if (bc_jhi[a_lev] != nullptr) { delete bc_jhi[a_lev]; }
  if (bc_klo[a_lev] != nullptr) { delete bc_klo[a_lev]; }
  if (bc_khi[a_lev] != nullptr) { delete bc_khi[a_lev]; }

  // Define and allocate the integer MultiFab that is the outside adjacent
  // cells of the problem domain.
  Box domainx(a_geom.Domain());
  domainx.grow(1,nghost);
  domainx.grow(2,nghost);
  Box box_ilo = amrex::adjCellLo(domainx,0,1);
  Box box_ihi = amrex::adjCellHi(domainx,0,1);

  Box domainy(a_geom.Domain());
  domainy.grow(0,nghost);
  domainy.grow(2,nghost);
  Box box_jlo = amrex::adjCellLo(domainy,1,1);
  Box box_jhi = amrex::adjCellHi(domainy,1,1);

  Box domainz(a_geom.Domain());
  domainz.grow(0,nghost);
  domainz.grow(1,nghost);
  Box box_klo = amrex::adjCellLo(domainz,2,1);
  Box box_khi = amrex::adjCellHi(domainz,2,1);

  // Note that each of these is a single IArrayBox so every
  // process has a copy of all BCs
  bc_ilo[a_lev] = new IArrayBox(box_ilo,2);
  bc_ihi[a_lev] = new IArrayBox(box_ihi,2);
  bc_jlo[a_lev] = new IArrayBox(box_jlo,2);
  bc_jhi[a_lev] = new IArrayBox(box_jhi,2);
  bc_klo[a_lev] = new IArrayBox(box_klo,2);
  bc_khi[a_lev] = new IArrayBox(box_khi,2);
}


MFIXBoundaryConditions::
MFIXBoundaryConditions ( int const a_max_level,
                         Vector<Geometry> const& a_geom,
                         Vector<IntVect> a_ref_ratio,
                         BCList& a_bc_list,
                         MFIXEmbeddedBoundaries& a_eb)
  : m_max_level(a_max_level), m_geom(a_geom), m_ref_ratio(a_ref_ratio)
  , m_bc_list(a_bc_list) , m_embedded_boundaries(a_eb)
{
  m_vel_bc_types["Dirichlet"] = {BCList::minf};
  m_vel_bc_types["Neumann"] = {BCList::pinf, BCList::pout};

  m_rho_bc_types["Dirichlet"] = {BCList::minf};
  m_rho_bc_types["Neumann"] = {BCList::pinf, BCList::pout};

  m_T_bc_types["Dirichlet"] = {BCList::minf, BCList::pinf};
  m_T_bc_types["Neumann"] = {BCList::pout};

  m_trac_bc_types["Dirichlet"] = {BCList::minf};
  m_trac_bc_types["Neumann"] = {BCList::pinf, BCList::pout};

  m_Xk_bc_types["Dirichlet"] = {BCList::minf, BCList::pinf};
  m_Xk_bc_types["Neumann"] = {BCList::pout};

}


void
MFIXBoundaryConditions::
Initialize ( const MFIXRegions& regions,
             MFIXFluidPhase& fluid,
             MFIXSolidsPhase& solids,
             MFIXDEM& dem,
             MFIXPIC& pic)
{

  // Read in verbosity flag
  { ParmParse pp;
    pp.query("mfix.verbose", m_verbose);
  }

  // Set flag to keep particles from leaving unless periodic.
  for (int dir(0); dir<3; ++dir) {
    if (m_geom[0].isPeriodic(dir)) {

      m_domain_bc[2*dir  ] = 0;
      m_domain_bc[2*dir+1] = 0;

    } else {
      m_domain_bc[2*dir  ] = 1;
      m_domain_bc[2*dir+1] = 1;
    }
  }

  // Default all sides of the domain to Neumann
  for (int dir=0; dir < 3; ++dir) {
    if (m_geom[0].isPeriodic(dir)) {

      m_ppe_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_ppe_hibc[dir] = amrex::LinOpBCType::Periodic;

      m_diff_vel_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_diff_vel_hibc[dir] = amrex::LinOpBCType::Periodic;

      m_diff_scal_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_diff_scal_hibc[dir] = amrex::LinOpBCType::Periodic;

      m_diff_temperature_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_diff_temperature_hibc[dir] = amrex::LinOpBCType::Periodic;

      m_diff_species_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_diff_species_hibc[dir] = amrex::LinOpBCType::Periodic;

      m_mac_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_mac_hibc[dir] = amrex::LinOpBCType::Periodic;

    } else {

      m_ppe_lobc[dir] = amrex::LinOpBCType::Neumann;
      m_ppe_hibc[dir] = amrex::LinOpBCType::Neumann;

      m_diff_vel_lobc[dir] = amrex::LinOpBCType::Dirichlet;
      m_diff_vel_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      m_diff_scal_lobc[dir] = amrex::LinOpBCType::Dirichlet;
      m_diff_scal_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      m_diff_temperature_lobc[dir] = amrex::LinOpBCType::Dirichlet;
      m_diff_temperature_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      m_diff_species_lobc[dir] = amrex::LinOpBCType::Dirichlet;
      m_diff_species_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      m_mac_lobc[dir] = amrex::LinOpBCType::Neumann;
      m_mac_hibc[dir] = amrex::LinOpBCType::Neumann;

    }
  }

  // Initialize the periodic pressure drop to zero in all directions
  for (int dir=0; dir < 3; ++dir) { m_delp[dir] = 0.0; }

  amrex::ParmParse pp_bc("bc");

  amrex::Real delp_in(0.0);
  pp_bc.query("delp", delp_in);

  if (delp_in != 0.0) {

    pp_bc.query("delp_dir", m_delp_dir);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0 <= m_delp_dir && m_delp_dir <= 2,
       "Direction of pressure drop (bc.delp_dir) must be specified.");

    m_delp[m_delp_dir] = delp_in;
  }

  constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

  // Flag to control if particles see a wall at pressure outflows
  int po_noParOut = 0; // default behavior for PO's -- letting particles exit the domain
  pp_bc.query("po_no_par_out", po_noParOut);

  // Get the list of region names used to define BCs
  std::vector<std::string> input_regions;
  pp_bc.queryarr("regions", input_regions);

  // Loop over BCs
  for(size_t bcv=0; bcv < input_regions.size(); bcv++) {

    amrex::Real volfrac_total(0.0);

    BC_t new_bc(&fluid);

    // Set the region for the initial condition.
    new_bc.region = regions.getRegion(input_regions[bcv]);

    if (new_bc.region == nullptr) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid bc region! " << input_regions[bcv];
    }

    // Get the BC type (MI/PI/PO/EB)
    std::string bc_type;
    pp_bc.get(input_regions[bcv], bc_type);

    // Convert the input string into the integers
    if( bc_type == "mi") {
      new_bc.type = BCList::minf;
    } else if( bc_type == "pi") {
      new_bc.type = BCList::pinf;
    } else if( bc_type == "po") {
      new_bc.type = BCList::pout;
    } else if (bc_type == "eb") {
      new_bc.type = BCList::eb;
    } else {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Unknown BC type: " << bc_type;
    }

    const auto plo = m_geom[0].ProbLoArray();
    const auto phi = m_geom[0].ProbHiArray();

    amrex::RealArray normal = {0.0, 0.0, 0.0};
    amrex::RealArray point;
    point[0] = new_bc.region->lo(0) + 0.5*(new_bc.region->hi(0)-new_bc.region->lo(0));
    point[1] = new_bc.region->lo(1) + 0.5*(new_bc.region->hi(1)-new_bc.region->lo(1));
    point[2] = new_bc.region->lo(2) + 0.5*(new_bc.region->hi(2)-new_bc.region->lo(2));

    std::string bc_region_prefix = "bc."+input_regions[bcv];
    amrex::ParmParse ppRegion(bc_region_prefix);

    int face = -1;

    if ( new_bc.type != BCList::eb ) {

      // This covers mass inflows (mi), pressure inflows (pi), and
      // pressure outflows (po). These all must align with a domain
      // extent.

      int sum_same_loc(0);
      for (int dir(0); dir<3; ++dir){

        int same_loc = amrex::Math::abs(new_bc.region->lo(dir) - new_bc.region->hi(dir)) < tolerance ? 1 : 0;

        if (same_loc) {

          if (amrex::Math::abs(new_bc.region->lo(dir) - plo[dir] ) < tolerance) {

            point[dir] = plo[dir]+1.0e-15;
            normal[dir] =  1.0;

            face = 2*dir;

            if (new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {
              m_ppe_lobc[dir] = amrex::LinOpBCType::Dirichlet;
              m_mac_lobc[dir] = amrex::LinOpBCType::Dirichlet;
              m_diff_vel_lobc[dir] = amrex::LinOpBCType::Neumann;
              m_diff_scal_lobc[dir] = amrex::LinOpBCType::Neumann;
            }

            if (new_bc.type == BCList::pout) {
              m_diff_temperature_lobc[dir] = amrex::LinOpBCType::Neumann;
              m_diff_species_lobc[dir] = amrex::LinOpBCType::Neumann;
            }

            if (new_bc.type == BCList::minf) {
              m_ppe_lobc[dir] = amrex::LinOpBCType::inflow;
            }


          } else if (amrex::Math::abs(new_bc.region->hi(dir) - phi[dir]) < tolerance) {

            point[dir] = phi[dir]-1.0e-15;
            normal[dir] = -1.0;

            face = 2*dir+1;

            if (new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {
              m_ppe_hibc[dir] = amrex::LinOpBCType::Dirichlet;
              m_mac_hibc[dir] = amrex::LinOpBCType::Dirichlet;
              m_diff_vel_hibc[dir] = amrex::LinOpBCType::Neumann;
              m_diff_scal_hibc[dir] = amrex::LinOpBCType::Neumann;
            }

            if (new_bc.type == BCList::pout) {
              m_diff_temperature_hibc[dir] = amrex::LinOpBCType::Neumann;
              m_diff_species_hibc[dir] = amrex::LinOpBCType::Neumann;
            }

            if (new_bc.type == BCList::minf) {
              m_ppe_hibc[dir] = amrex::LinOpBCType::inflow;
            }

          } else {

            reporter::Log(reporter::Error,__FILE__, __LINE__) << " Error:"
              << "Flow BCs must be located on domain extents.\n"
              << "BC Name: " << input_regions[bcv] << "  "
              << "Invalid direction: " << dir;
          }

          // Flag that the level set should see this domain extent as walls
          // for particle collisions.
          if ( new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
              (new_bc.type == BCList::pout && po_noParOut == 1)) {
            m_ls_wall_flag.set(face);
          }
        }

        sum_same_loc += same_loc;
      }

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( sum_same_loc == 1, "Invalid bc region!");
    }

    // Store the BC ID for quick look-up when setting BC types
    if (face == 0) { m_bc_lo[0].push_back(bcv); }
    else if (face == 1) { m_bc_hi[0].push_back(bcv); }
    else if (face == 2) { m_bc_lo[1].push_back(bcv); }
    else if (face == 3) { m_bc_hi[1].push_back(bcv); }
    else if (face == 4) { m_bc_lo[2].push_back(bcv); }
    else if (face == 5) { m_bc_hi[2].push_back(bcv); }

    m_dir.push_back(face);

    // Enforce the boundary for pressure outflows if specified.
    if (new_bc.type == BCList::pout && (po_noParOut == 0))
    { m_domain_bc[face] = 0; }

    // Get EB data.
    if (new_bc.type == BCList::eb) {

      std::string bc_region_eb_prefix = "bc."+input_regions[bcv]+".eb";
      ParmParse ppEB(bc_region_eb_prefix);

      std::string eb_bc_T = "adiabatic";
      if (ppEB.queryAdd("temperature", eb_bc_T)) {

        if ( amrex::toLower(eb_bc_T).compare("adibatic") == 0) {

          // Nothing to do. EB walls are adibatic.

          // This is for debugging and could be removed.
          reporter::Log(reporter::Status)
            << "EB is adibatic:  bc." << input_regions[bcv];

        } else {

          // Set flag for non-adibatic EB walls
          eb_parms().set_has_temperature();

          // Constant temperature walls
          if (amrex::toLower(eb_bc_T).compare("constant") == 0) {

            if (new_bc.eb.m_temperature.get(ppEB, "temperature.constant")) {

              reporter::Log(reporter::Error,__FILE__, __LINE__)
                << "Failed to process boundary condition EB temperature!\n"
                << "BC region: " << input_regions[bcv] << '\n'
                << bc_region_eb_prefix << ".temperature.constant";
            }

            new_bc.eb.set_constant_temperature();

            // This is for debugging and could be removed.
            reporter::Log(reporter::Status)
              << "EB with constant temperature:  bc." << input_regions[bcv]
              << "\n  Temperature of wall .................... " << new_bc.eb.temperature(0.);

          // 1D Conjugate heat transfer walls
          } else if (amrex::toLower(eb_bc_T).compare("cht-1d") == 0) {

            if (!ppEB.queryAdd("temperature.CHT-1D.T_env", new_bc.eb.T_env)) {

              // Log that the default is being used.
              reporter::Log(reporter::Status)
                << "EB conjugate heat transfer 1D input not specified. Setting default:\n"
                << "  " << bc_region_eb_prefix << ".temperature.CHT-1D.T_env" << new_bc.eb.T_env;
            }
            if (!ppEB.queryAdd("temperature.CHT-1D.h_int", new_bc.eb.h_int)) {

              // Log that the default is being used.
              reporter::Log(reporter::Status)
                << "EB conjugate heat transfer 1D input not specified. Setting default:\n"
                << "  " << bc_region_eb_prefix << ".temperature.CHT-1D.h_int" << new_bc.eb.h_int;
            }
            if (!ppEB.queryAdd("temperature.CHT-1D.h_ext", new_bc.eb.h_ext)) {

              // Log that the default is being used.
              reporter::Log(reporter::Status)
                << "EB conjugate heat transfer 1D input not specified. Setting default:\n"
                << "  " << bc_region_eb_prefix << ".temperature.CHT-1D.h_ext" << new_bc.eb.h_ext;
            }
            if (!ppEB.queryAdd("temperature.CHT-1D.wall.conductivity", new_bc.eb.k_wall)) {

              // Log that the default is being used.
              reporter::Log(reporter::Status)
                << "EB conjugate heat transfer 1D input not specified. Setting default:\n"
                << "  " << bc_region_eb_prefix << ".temperature.CHT-1D.wall.conductivity" << new_bc.eb.k_wall;
            }
            if (!ppEB.queryAdd("temperature.CHT-1D.wall.thickness", new_bc.eb.L_wall)) {

              // Log that the default is being used.
              reporter::Log(reporter::Status)
                << "EB conjugate heat transfer 1D input not specified. Setting default:\n"
                << "  " << bc_region_eb_prefix << ".temperature.CHT-1D.wall.thickness" << new_bc.eb.L_wall;
            }

            ppEB.queryAdd("temperature.CHT-1D.phase_averaged_temperature",
                new_bc.eb.phase_averaged_temperature);

            if (!dem.solve() && !pic.solve() && new_bc.eb.phase_averaged_temperature ) {

              // Log that the default is being used.
              reporter::Log(reporter::Warning)
                << "EB conjugate heat transfer 1D cannot use phase averaged temperature\n"
                << "when no solid model is enabled!\n"
                << "Set override: " << bc_region_eb_prefix
                << ".temperature.CHT-1D.phase_averaged_temperature = false";

                new_bc.eb.phase_averaged_temperature = false;
            }

            // Set flag indicating we need to compute the averaged particle temperature
            if ( new_bc.eb.phase_averaged_temperature ) {
              eb_parms().set_needs_avg_particle_temperature();
            }

            new_bc.eb.set_conjugate_heat_transfer_1D();

            // This is for debugging and could be removed.
            reporter::Log(reporter::Status)
              << "EB conjugate heat transfer 1D settings:  bc." << input_regions[bcv]
              << "\n  Temperature of environment ............. " << new_bc.eb.T_env
              << "\n  Internal heat transfer coefficient ..... " << new_bc.eb.h_int
              << "\n  External heat transfer coefficient ..... " << new_bc.eb.h_ext
              << "\n  Wall conductivity ...................... " << new_bc.eb.k_wall
              << "\n  Wall thickness ......................... " << new_bc.eb.L_wall
              << "\n  Phase averaged temperature ............. " << new_bc.eb.phase_averaged_temperature;

          } else {

            // Flag error -- Unknown EB wall temperature type
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Unknown EB wall temperature type: " << eb_bc_T;
          }
        } // end non-adibatic wall options
      } // end EB settings for temperature

      if (ppEB.contains("normal")) {
        new_bc.eb.has_normal = 1;

        ppEB.getarr("normal", new_bc.eb.normal, 0, 3);
        const amrex::Real nx = new_bc.eb.normal[0];
        const amrex::Real ny = new_bc.eb.normal[1];
        const amrex::Real nz = new_bc.eb.normal[2];

        const amrex::Real nmag = std::sqrt(nx*nx+ny*ny+nz*nz);
        AMREX_ALWAYS_ASSERT(nmag > 0.);

        new_bc.eb.normal[0] = nx / nmag;
        new_bc.eb.normal[1] = ny / nmag;
        new_bc.eb.normal[2] = nz / nmag;

        amrex::Real tol_deg(0.);
        ppEB.query("normal_tol", tol_deg);

        new_bc.eb.normal_tol = tol_deg*M_PI/amrex::Real(180.);
      }
    }

    // Get fluid data.
    if (fluid.solve()) {

      std::string bc_region_fluid_prefix = "bc."+input_regions[bcv]+"."+fluid.names(0);
      amrex::ParmParse ppFluid(bc_region_fluid_prefix);

      // Mass inflows need fluid velocity and volume fraction.
      if (new_bc.type == BCList::minf || new_bc.type == BCList::eb) {

        int has_flow(0);

        if (ppFluid.contains("velocity")) {

           if (new_bc.fluid.velocity.get(ppFluid, "velocity")) {

             reporter::Log(reporter::Error,__FILE__, __LINE__)
               << "Failed to process boundary condition fluid velocity.\n"
               << "BC region: " << input_regions[bcv] << '\n'
               << "Please correct the input deck.";
           }

           has_flow = 1;

        } else if (ppFluid.contains("volflow")) {

           if (new_bc.fluid.volflow.get(ppFluid, "volflow")) {

             reporter::Log(reporter::Error,__FILE__, __LINE__)
               << "Failed to process boundary condition fluid volflow.\n"
               << "BC region: " << input_regions[bcv] << '\n'
               << "Please correct the input deck.";
           }

           has_flow = 1;

        } else if (ppFluid.contains("massflow")) {

          if (new_bc.fluid.massflow.get(ppFluid, "massflow")) {

            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Failed to process boundary condition fluid velocity.\n"
              << "BC region: " << input_regions[bcv] << '\n'
              << "Please correct the input deck.";
          }

           has_flow = 1;
        }

        if (has_flow) {

          ppFluid.get("volfrac", new_bc.fluid.volfrac);
          volfrac_total += new_bc.fluid.volfrac;

          if (new_bc.type == BCList::eb) {
            new_bc.fluid.flow_thru_eb = 1;
            eb_parms().set_has_flow();
          }
        }
      } // MI or EB

      // Read in density, temperature and species BC inputs only if BC type is
      // minf, pinf or flow_through_eb
      if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
          (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb)) {

        // Read in fluid density
        if (new_bc.fluid.density.get(ppFluid, "density")) {

          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Failed to process boundary condition fluid density.\n"
            << "BC region: " << input_regions[bcv] << '\n'
            << "Please correct the input deck.";
        }


        if (new_bc.fluid.temperature.get(ppFluid, "temperature")) {

          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Failed to process boundary condition fluid temperature.\n"
            << "BC region: " << input_regions[bcv] << '\n'
            << "Please correct the input deck.";
        }

        // Get species data.
        if (fluid.solve_species()) {

          int const nspecies_g = fluid.nspecies();
          new_bc.fluid.species.resize(nspecies_g, FVAR_({1}));

          std::string bc_region_fluid_species_prefix = bc_region_fluid_prefix+".species";
          ParmParse ppSpecies(bc_region_fluid_species_prefix);

          // Auxiliary variable to check that species sum up to 1
          Real massfrac_total(0);

          for (int n(0); n < nspecies_g; n++) {

            std::string species_name = fluid.species_names(n);
            if (new_bc.fluid.species[n].get(ppSpecies, species_name)) {
              reporter::Log(reporter::Error,__FILE__, __LINE__)
                << "Failed to process boundary condition fluid species.\n"
                << "BC region: " << input_regions[bcv] << "   "
                << "species: " << species_name;
            }

            massfrac_total += new_bc.fluid.get_species(n);
          }

          if ( !almostEqual(massfrac_total,1.0) ) {
            reporter::Log(reporter::Warning)
              << "Fluid mass fractions do not sum to one.\n"
              << "BC region: " << input_regions[bcv]
              << "  Total volume fraction: " << massfrac_total;
          }
        }
      }

      // Read in fluid pressure BC inputs only if BC type is pinf or pout
      if (new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {

        // Read in fluid pressure
        new_bc.fluid.pressure_defined = ppFluid.query("pressure", new_bc.fluid.pressure);

        if (!new_bc.fluid.pressure_defined) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Pressure BCs must have pressure defined.\n"
            << "BC region: " << input_regions[bcv] << '\n'
            << "Please correct the input deck.";
        }
      }

      // ********************************************************************
      // Check BCs consistency depending on fluid constraint type
      // ********************************************************************
      if ( fluid.constraint.isIncompressibleFluid() ) {

        if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
            (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb)) {

          if (!new_bc.fluid.density.is_defined()) {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Fluid density is required for IncompressibleFLuid constraint.\n"
              << "BC region: " << input_regions[bcv] << "   "
              << "Please correct the input deck.";
          }
        }

      } else if ( fluid.constraint.isIdealGasOpenSystem() ) {

        if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
            (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb)) {

          if (new_bc.fluid.density.is_defined()) {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Fluid density should not be defined with IdealGasOpenSystem constraint.\n"
              << "BC region: " << input_regions[bcv] << '\n'
              << "Please correct the input deck.";
          }

          if (!new_bc.fluid.temperature.is_defined()) {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Fluid temperature is required for IdealGasOpenSystem constraint.\n"
              << "BC region: " << input_regions[bcv] << '\n'
              << "Please correct the input deck.";
          }
        }

      } else if ( fluid.constraint.isIdealGasClosedSystem() ) {

        if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Open boundaries are not available with IdealGasClosedSystem constraint.\n"
              << "BC region: " << input_regions[bcv] << '\n'
              << "Please correct the input deck.";
        }

        if (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb) {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Flow through EB is not available with IdealGasClosedSystem constraint.\n"
              << "BC region: " << input_regions[bcv] << '\n'
              << "Please correct the input deck.";
        }
      }

    }

    // Get solid data.
    if(dem.solve() || pic.solve()) {

      // Get the list of solids used in defining the BC region
      std::vector<std::string> solids_names;
      ppRegion.queryarr("solids", solids_names);

      for(size_t lcs(0); lcs < solids_names.size(); ++ lcs){

        SOLIDS_t new_solid;

        new_solid.name = solids_names[lcs];
        new_solid.phase = solids.name_to_phase(solids_names[lcs]);

        std::string bc_region_solid_prefix = "bc."+input_regions[bcv]+"."+solids_names[lcs];
        ParmParse ppSolid(bc_region_solid_prefix);

        ppSolid.get("volfrac", new_solid.volfrac);
        volfrac_total += new_solid.volfrac;

        bool const is_mi(new_bc.type == BCList::minf);
        bool const is_eb(new_bc.type == BCList::eb);

        if (new_solid.volfrac > tolerance) {

          if (is_mi && dem.solve()) {

            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Particle inflow boundaries are not supported for DEM solids.\n"
              << "BC region: " << input_regions[bcv] << '\n'
              << "Please correct the input deck.";

          } else if (is_eb || (is_mi && pic.solve()) ) {

            int const dir_lohi(face%2); // lo=0, hi=1
            int const dir((face-dir_lohi)/((int)2)); // dir=0,1,2 (x,y,z)

            if (ppSolid.contains("velocity")) {

              amrex::Vector<amrex::Real> vel_in;
              ppSolid.queryarr("velocity", vel_in);

              const int rcomps = vel_in.size();

              if (rcomps == 3) {
                new_solid.velocity = vel_in;
                new_solid.velmag = std::sqrt(vel_in[0]*vel_in[0]
                                           + vel_in[1]*vel_in[1]
                                           + vel_in[2]*vel_in[2]);

              } else if (rcomps == 1) {
                new_solid.velmag = amrex::Math::abs(vel_in[0]);

                new_solid.velocity.push_back( (dir == 0 ? vel_in[0] : 0.) );
                new_solid.velocity.push_back( (dir == 1 ? vel_in[0] : 0.) );
                new_solid.velocity.push_back( (dir == 2 ? vel_in[0] : 0.) );

              } else {
                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "Invalid solids velocity specified.\n"
                  << "BC region: " << input_regions[bcv] << "   "
                  << "solids: " << solids_names[lcs];
              }

            } else if (ppSolid.query("volflow", new_solid.volflow)) {

              // +1 low side; -1 high side
              Real const r_normal_dir = (dir_lohi == 0) ? 1. : -1.;

              new_solid.velocity.push_back( ( dir == 0 ? r_normal_dir : 0.) );
              new_solid.velocity.push_back( ( dir == 1 ? r_normal_dir : 0.) );
              new_solid.velocity.push_back( ( dir == 2 ? r_normal_dir : 0.) );

            } else {
                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "Solids velocity not specified.\n"
                  << "BC region: " << input_regions[bcv] << "   "
                  << "solids: " << solids_names[lcs];
            }

            // This probably needs a check if we are solving energy equations
            // and if so, require temperature be defined.
            ppSolid.query("temperature", new_solid.temperature);

            { // Get information about diameter distribution.
              if( new_solid.diameter.define(bc_region_solid_prefix, "diameter",
                    new_solid.get_diameter_data(), dem.solve()) ) {
                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "BC region " + input_regions[bcv] + " has invalid diameter parameters "
                  << "for solids " + solids_names[lcs];
              }
            }

            if (solids.solve_species()) {

              const int nspecies_s = solids.nspecies();
              new_solid.species.resize(nspecies_s);

              std::string bc_region_solid_species_prefix = bc_region_solid_prefix+".species";
              ParmParse ppSpecies(bc_region_solid_species_prefix);

              Real total_mass_fraction(0);

              for (int n(0); n < solids.nspecies(); n++) {

                new_solid.species[n].mass_fraction = 0.;

                std::string species_name = solids.species_names(n);
                ppSpecies.query(species_name, new_solid.species[n].mass_fraction);

                total_mass_fraction += new_solid.species[n].mass_fraction;
              }

              // Sanity check that the input species mass fractions sum up to 1
              if (!almostEqual(total_mass_fraction,1.) ) {

                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "Solids mass fractions do not sum to 1.\n"
                  << "BC region: " << input_regions[bcv] << "   "
                  << "solids: " << solids_names[lcs];
              }
            } // end solve species

            { // Get information about density distribution
              if (solids.density_model() == MFIXSolidsPhase::DensityModel::SpeciesMixture)
              {
                amrex::Real p_density_inv(0.);
                for (int n(0); n < solids.nspecies(); n++) {
                  p_density_inv += new_solid.species[n].mass_fraction / solids.species_density(n);
                }
                amrex::Real p_density = 1./p_density_inv;

                std::string density_type("");
                const bool contains_density = ppSolid.query("density", density_type);

                if (contains_density) {
                  if (amrex::toLower(density_type) == "constant") {
                    amrex::Real input_density(0.);
                    ppSolid.get("density.constant", input_density);

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
                        << "BC region: " << input_regions[bcv] << "   "
                        << "solids: " << solids_names[lcs];
                    }
                  } else {
                    reporter::Log(reporter::Error,__FILE__, __LINE__)
                      << "Invalid particle density parameters.\n"
                      << "When solids density model is a mixture of species densities\n"
                      << "the BC solids density should either not be specified"
                      << " or specified as constant.\n"
                      << "BC region: " << input_regions[bcv] << "   "
                      << "solids: " << solids_names[lcs];
                  }
                } else {
                  ppSolid.add("density", std::string("constant"));
                  ppSolid.add("density.constant", amrex::Real(p_density));
                }
              }

              if( new_solid.density.define(bc_region_solid_prefix, "density") ) {
                reporter::Log(reporter::Error,__FILE__, __LINE__)
                  << "BC region " + input_regions[bcv] + " has invalid density parameters "
                  << "for solids " + solids_names[lcs];
              }

            }

          } // End EB inflow
        } // volfrac > tolerance

        new_bc.solids.push_back(new_solid);
      }
    }

    if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
         (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb)) {

      if (fluid.solve() && (dem.solve() || pic.solve())) {

        if ( !almostEqual(volfrac_total,1.0) ) {
              reporter::Log(reporter::Warning)
                << "Fluid and particle volume fractions do not sum to one.\n"
                << "BC region: " << input_regions[bcv]
                << "  Total volume fraction: " << volfrac_total;
        }
      }
    }

    m_bc.push_back(new_bc);
  }

  if (fluid.solve() &&
      (fluid.constraint.isIncompressibleFluid() ||
       fluid.constraint.isIdealGasOpenSystem() )) {

    Vector<Real> bc_values(0);

    for (int bcv(0); bcv < m_bc.size(); ++bcv) {

      if (m_bc[bcv].type == BCList::pout)
        bc_values.push_back(m_bc[bcv].fluid.pressure);
    }

    if (bc_values.size() > 0) {
      const Real max_val = *std::max_element(bc_values.begin(), bc_values.end());
      const Real min_val = *std::min_element(bc_values.begin(), bc_values.end());

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(max_val-min_val) < 1.e-15,
        "Error: pressure values at pout boundaries differ");

      if (fluid.constraint.isIdealGasOpenSystem() ) {

        if (fluid.p_therm_defined()) {
          const Real p_therm = fluid.thermodynamic_pressure();

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(p_therm-min_val) < 1.e-15 &&
              std::abs(p_therm-max_val) < 1.e-15 && p_therm > 1.e-15,
              "Error: either pout input pressures differ from fluid thermodynamic "
              "pressure, or given thermodynamic pressure is zero or negative");

        } else {

          fluid.set_thermodynamic_pressure(max_val);
          fluid.set_p_therm_defined(1);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.thermodynamic_pressure() > 1.e-15,
              "Error: fluid thermodynamic pressure is zero or negative");

        }
      }
    }
  }

  // This completes the copyAsync to Gpu::DeviceVectors and
  // sets pointers for runtime use.
  for (int bcv(0); bcv<m_bc.size(); ++bcv) {
    for (int m(0); m<m_bc[bcv].solids.size(); ++m) {
      m_bc[bcv].solids[m].copyAsync();
    }
  }

} // END Initialize


amrex::Real& MFIXBoundaryConditions::
get_bc_area (const int a_bcv)
{
  AMREX_ASSERT( a_bcv < m_bc.size() );
  AMREX_ASSERT( m_bc[a_bcv].type == BCList::pinf ||
                m_bc[a_bcv].type == BCList::pout ||
                m_bc[a_bcv].type == BCList::minf ||
                m_bc[a_bcv].type == BCList::eb );

  return m_area[a_bcv];
}

amrex::Real& MFIXBoundaryConditions::
get_bc_fab_area (const int a_bcv, std::pair<int,int> a_index)
{
  AMREX_ASSERT( a_bcv < m_bc.size() );
  AMREX_ASSERT( m_bc[a_bcv].type == BCList::pinf ||
                m_bc[a_bcv].type == BCList::pout ||
                m_bc[a_bcv].type == BCList::minf ||
                m_bc[a_bcv].type == BCList::eb );

  return m_fab_area[a_bcv][a_index];
}
