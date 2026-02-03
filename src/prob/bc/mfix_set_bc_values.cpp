#include <mfix_run_on.H>
#include <mfix_reporter.H>
#include <mfix_bc.H>
#include <mfix_fluid.H>

using namespace amrex;

void
MFIXBoundaryConditions::
set_velocity_bc_values (Real time_in)
{

  m_h_bc_u_g.resize(m_bc.size());
  m_h_bc_v_g.resize(m_bc.size());
  m_h_bc_w_g.resize(m_bc.size());

  m_bc_u_g.resize(m_bc.size());
  m_bc_v_g.resize(m_bc.size());
  m_bc_w_g.resize(m_bc.size());


  for (unsigned bcv(0); bcv < m_bc.size(); ++bcv) {

    const BC_t& bc = m_bc[bcv];

    if (bc.type == BCList::minf) {

      if (bc.fluid.velocity.size() == 3) {

          RealVect bc_vels = bc.fluid.get_velocity(time_in);

          m_h_bc_u_g[bcv] = bc_vels[0];
          m_h_bc_v_g[bcv] = bc_vels[1];
          m_h_bc_w_g[bcv] = bc_vels[2];

      } else {

        Real vel_mag(0.);

        if (bc.fluid.velocity.is_defined()) {

          AMREX_ASSERT( bc.fluid.velocity.size() == 1 );
          vel_mag = bc.fluid.get_velocity_mag(time_in);

        } else {

          std::vector<std::string> input_regions;
          amrex::ParmParse pp("bc");
          pp.queryarr("regions", input_regions);

          Real volflow = std::numeric_limits<Real>::max();

          // massflow --> volflow
          if (bc.fluid.massflow.is_defined()) {

            Real massflow = bc.fluid.get_massflow(time_in);
            Real density = bc.fluid.get_density(time_in);

            if (amrex::almostEqual(density, 0.)) {

              // Error check
              reporter::Log(reporter::Error,__FILE__, __LINE__)
                << "Unable to convert BC mass flow to volumetric flow!\n"
                << "BC region: " << input_regions[bcv] << "  "
                << "density: " << density << '\n'
                << "Please correct the input deck.";
            }

            volflow = massflow / density;

          } else if (bc.fluid.volflow.is_defined()) {

            volflow = bc.fluid.get_volflow(time_in);

          } else {

            // Error check -- HOW DID WE GET HERE?
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Invalid boundary flow!\n"
              << "BC region: " << input_regions[bcv] << '\n'
              << "Please correct the input deck.";
          }

          Real area = get_bc_area(bcv);
          area *= bc.fluid.volfrac;

          if (amrex::almostEqual(area, 0.)) {

            // Error check
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Unable to convert BC volumetric to velocity!\n"
              << "BC region: " << input_regions[bcv] << "  "
              << "BC area: " << area << '\n'
              << "Please correct the input deck.";

          }
          vel_mag = volflow / area;
        }

        // 0:x-lo, 1:x-hi, 2:y-lo, 3:y-hi, 4:z-lo, 5:z-hi
        int const face = get_dir(bcv);

        int const dir_lohi(face%2); // lo=0, hi=1
        int const dir((face-dir_lohi)/((int)2)); // dir=0,1,2 (x,y,z)

        vel_mag *= ((dir_lohi == 0) ? 1. : -1.);

        m_h_bc_u_g[bcv] = ((dir == 0) ? vel_mag : 0.);
        m_h_bc_v_g[bcv] = ((dir == 1) ? vel_mag : 0.);
        m_h_bc_w_g[bcv] = ((dir == 2) ? vel_mag : 0.);

      }

    } else {

      m_h_bc_u_g[bcv] = 1e50;
      m_h_bc_v_g[bcv] = 1e50;
      m_h_bc_w_g[bcv] = 1e50;

    }

  }

  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_u_g.begin(), m_h_bc_u_g.end(), m_bc_u_g.begin());

  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_v_g.begin(), m_h_bc_v_g.end(), m_bc_v_g.begin());

  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_w_g.begin(), m_h_bc_w_g.end(), m_bc_w_g.begin());

  Gpu::synchronize();
}


void
MFIXBoundaryConditions::
set_energy_bc_values ( Real a_time,
                       const MFIXFluidPhase& a_fluid)
{

  int const has_temperature = a_fluid.solve_enthalpy() ||
      a_fluid.constraint.isIdealGas();

  if (has_temperature) {
    m_h_bc_T.resize(m_bc.size());
    m_bc_T.resize(m_bc.size());
  }

  int const has_enthalpy = a_fluid.solve_enthalpy();

  if (has_enthalpy) {
    m_h_bc_h.resize(m_bc.size());
    m_bc_h.resize(m_bc.size());
  }

  const auto fluid_props = a_fluid.props.data<RunOn::Host>();

  for(unsigned bcv(0); bcv < m_bc.size(); ++bcv) {

    const BC_t& bc = m_bc[bcv];

    if (bc.type == BCList::minf || bc.type == BCList::pinf ||
        (bc.type == BCList::eb && bc.fluid.flow_thru_eb)) {

      AMREX_ASSERT( bc.fluid.temperature.is_defined() );
      const Real Tf = bc.fluid.get_temperature(a_time);

      m_h_bc_T[bcv] = Tf;

      if (has_enthalpy) {

        if (!a_fluid.isMixture()) {

          m_h_bc_h[bcv] = fluid_props.enthalpy(m_h_bc_T[bcv], nullptr);

        } else {

          m_h_bc_h[bcv] = 0.0;

          for (int n(0); n < a_fluid.nspecies(); n++) {

            const Real X_gk = bc.fluid.get_species(n, a_time);
            m_h_bc_h[bcv] += X_gk*fluid_props.enthalpy(n,Tf);
          }
        }

      }

    // BC is not minf, pinf, eb
    } else {

      if (has_temperature) { m_h_bc_T[bcv] = 1e50; }
      if (has_enthalpy) { m_h_bc_h[bcv] = 1e50; }
    }
  }

  if (has_temperature) {
    Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_T.begin(), m_h_bc_T.end(), m_bc_T.begin());
  }

  if (has_enthalpy) {
    Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_h.begin(), m_h_bc_h.end(), m_bc_h.begin());
  }

  Gpu::synchronize();
}


void MFIXBoundaryConditions::
set_epf_values ( int const a_solve_fluid )
{
  m_h_bc_epf.resize(m_bc.size());
  m_bc_epf.resize(m_bc.size());

  for(unsigned bcv(0); bcv < m_bc.size(); ++bcv) {

    const BC_t& bc = m_bc[bcv];

    m_h_bc_epf[bcv] = a_solve_fluid ? bc.fluid.volfrac : 1e50;
  }

  Gpu::copy(Gpu::hostToDevice, m_h_bc_epf.begin(), m_h_bc_epf.end(), m_bc_epf.begin());
}


void MFIXBoundaryConditions::
set_pressure_values ( int const a_solve_fluid )
{
  m_h_bc_p.resize(m_bc.size());
  m_bc_p.resize(m_bc.size());

  for(unsigned bcv(0); bcv < m_bc.size(); ++bcv) {

    const BC_t& bc = m_bc[bcv];

    m_h_bc_p[bcv] = a_solve_fluid ? bc.fluid.pressure : 1e50;
  }

  Gpu::copy(Gpu::hostToDevice, m_h_bc_p.begin(), m_h_bc_p.end(), m_bc_p.begin());
}


void
MFIXBoundaryConditions::
set_density_bc_values ( Real time_in )
{
  m_h_bc_rho.resize(m_bc.size());
  m_bc_rho.resize(m_bc.size());

  for(unsigned bcv(0); bcv < m_bc.size(); ++bcv) {

    const BC_t& bc = m_bc[bcv];

    if (bc.type == BCList::minf || bc.type == BCList::pinf ||
        (bc.type == BCList::eb && bc.fluid.flow_thru_eb)) {

      m_h_bc_rho[bcv] = bc.fluid.get_density(time_in);

    } else {
      m_h_bc_rho[bcv] = 1e50;
    }
  }

  Gpu::copy(Gpu::hostToDevice, m_h_bc_rho.begin(), m_h_bc_rho.end(), m_bc_rho.begin());
}


void
MFIXBoundaryConditions::
set_tracer_bc_values (Real /*time_in*/,
                      const MFIXFluidPhase& fluid)
{
  m_h_bc_tracer.resize(m_bc.size());
  m_bc_tracer.resize(m_bc.size());

  // HACK -- BC tracer is constant given current implementation.
  // This was copied over from the mfix_set_tracer_bcs routine.
  const Real trac0 = fluid.tracer();

  for(unsigned bcv(0); bcv < m_bc.size(); ++bcv) {

    const BC_t& bc = m_bc[bcv];

    if (bc.type == BCList::minf || bc.type == BCList::pinf ||
        (bc.type == BCList::eb && bc.fluid.flow_thru_eb) ) {

      m_h_bc_tracer[bcv] = trac0;
    } else {
      m_h_bc_tracer[bcv] = 1e50;
    }
  }

  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_tracer.begin(), m_h_bc_tracer.end(), m_bc_tracer.begin());

  Gpu::synchronize();
}

void
MFIXBoundaryConditions::
set_species_bc_values ( Real time_in,
                        int const a_nspecies )
{
  int const stride(m_bc.size());

  m_h_bc_Xk.resize(stride*a_nspecies);
  m_bc_Xk.resize(stride*a_nspecies);

  for (int n(0); n < a_nspecies; n++) {

    for(int bcv(0); bcv < m_bc.size(); ++bcv) {

      const BC_t& bc = m_bc[bcv];

      if ( bc.type == BCList::minf || bc.type == BCList::pinf ||
         ( bc.type == BCList::eb && bc.fluid.flow_thru_eb) ) {

        AMREX_ASSERT( bc.fluid.species.size() == a_nspecies );

        m_h_bc_Xk[bcv + n*stride] = bc.fluid.get_species(n, time_in);


      } else {

        m_h_bc_Xk[bcv + n*stride] = 1e50;

      }
    } // bcv
  } // nspecies

  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_Xk.begin(),
      m_h_bc_Xk.end(), m_bc_Xk.begin());

  Gpu::synchronize();
}
