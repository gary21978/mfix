#include <limits>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <AMReX_AmrParGDB.H>

#include <mfix_pc_generator.H>

#include <mfix_pc_hex_close_pack_K.H>
#include <mfix_pc_n-cube_per_fill_K.H>
#include <mfix_pc_random_fill_dem_K.H>
#include <mfix_pc_random_fill_pic_K.H>

#include <mfix_reporter.H>
#include <mfix_solids.H>

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                !
//                                                                      !
//  Purpose: Generate particle configuration based on maximum particle  !
//           radius and filling from top to bottom within specified     !
//           bounds                                                     !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

using namespace amrex;


ParticlesGenerator::ParticlesGenerator ( const RealVect& a_plo,
                                         const RealVect& a_dx,
                                         amrex::RealBox const* const a_p_regions,
                                         int const a_regions_nb,
                                         int const a_allow_overlap,
                                         const RealBox& a_ic_region,
                                         const SOLIDS_t& a_ic_solid,
                                         std::string const a_ic_pack_type,
                                         int const a_has_granular_temperature,
                                         int const a_is_dem,
                                         int const a_is_pic,
                                         int const a_is_cg_dem,
                                         Real const a_multi_particle_volume)
  : m_plo(a_plo)
  , m_dx(a_dx)
  , m_ic_region(a_ic_region)
  , m_ic_solid(a_ic_solid)
  , m_is_dem(a_is_dem)
  , m_is_pic(a_is_pic)
  , m_is_cg_dem(a_is_cg_dem)
  , m_multi_particle_volume(a_multi_particle_volume)
  , m_is_hex_close_pack( a_ic_pack_type.compare("hcp") == 0 )
  , m_is_random_pack( a_ic_pack_type.compare("random") == 0 )
  , m_is_pseudo_random_pack( a_ic_pack_type.compare("pseudo_random") == 0 )
  , m_cube_base(-1)
  , m_has_granular_temperature(a_has_granular_temperature)
  , m_p_regions(a_p_regions)
  , m_regions_nb(a_regions_nb)
  , m_allow_overlap(a_allow_overlap)
{
  amrex::ignore_unused(m_is_pic);

  // Get the X-ncube integer input
  size_t const dash( a_ic_pack_type.find_first_of("-") );
  if ( dash != std::string::npos)
  { try { m_cube_base =  std::stoi(a_ic_pack_type.substr(0,dash)); }
    catch (...) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Unable to determine arrangement from cube packing: " << a_ic_pack_type << "\n"
        << "Input format: <int>-cube\n"
        << "Please correct the input deck.";
    }
  }
  else if (a_ic_pack_type.compare("oneper") == 0) { m_cube_base = 1; }
  else if (a_ic_pack_type.compare("eightper") == 0) { m_cube_base = 2; }
}



int
ParticlesGenerator::generate ( Box const& a_bx,
                               int const a_id, int const a_cpu,
                               ParticleTileType& particles)
{
  int num_particles(0);

  // Call generate with specific positions generator
  if (m_is_dem && m_is_hex_close_pack) {

    Hex_ClosePack hex_close_pack(m_plo, m_dx);

    hex_close_pack.setup(a_bx, m_ic_region,
      m_ic_solid.diameter.get_mean(), m_ic_solid.volfrac);

    num_particles = generate(a_id, a_cpu, particles, hex_close_pack);

  } else if ( m_is_random_pack || m_is_pseudo_random_pack ) {

    if (m_is_dem) {

      amrex::Gpu::HostVector<amrex::Real> h_pos, h_rad;
      amrex::Gpu::DeviceVector<amrex::Real> d_pos, d_rad;

      RandomFill_DEM random_fill_dem(m_plo, m_dx,
              *m_ic_solid.get_diameter(), m_is_random_pack);

      random_fill_dem.setup(a_bx, m_ic_solid, m_multi_particle_volume,
              h_pos, h_rad, d_pos, d_rad);

      num_particles = generate(a_id, a_cpu, particles, random_fill_dem);

    } else { AMREX_ASSERT(m_is_pic);

      RandomFill_PIC random_fill_pic(m_plo, m_dx, *m_ic_solid.get_diameter(),
        m_is_random_pack);

      random_fill_pic.setup(a_bx, m_ic_solid, m_multi_particle_volume);

      num_particles = generate(a_id, a_cpu, particles, random_fill_pic);
    }

  } else if (m_cube_base > 0) {

    nCubePer_Fill n_cube_per_fill(m_plo, m_dx, *m_ic_solid.get_diameter(), m_cube_base);
    n_cube_per_fill.setup(a_bx);

    num_particles = generate(a_id, a_cpu, particles, n_cube_per_fill);

  } else {

    amrex::Abort("Unknown particle generator fill type");
  }
  return num_particles;
}


template <typename F1>
int ParticlesGenerator::generate ( int const a_id, int const a_cpu,
                                   ParticleTileType& particles,
                                   F1 a_generator)
{
  int particles_count = a_generator.get_particles_number();

  const int current_size = particles.numParticles();

  particles.resize(current_size + particles_count);


  {

    // Setup particle diameters parameters
    const auto& diameter = m_ic_solid.diameter;

    // Setup particle densities parameters
    const auto& density = m_ic_solid.density;

    // Get particles initial velocity
    const Real ic_u_s = m_ic_solid.velocity[0];
    const Real ic_v_s = m_ic_solid.velocity[1];
    const Real ic_w_s = m_ic_solid.velocity[2];

    const int has_granular_temperature = m_has_granular_temperature;

    const Real mpv = m_multi_particle_volume;

    auto& aos = particles.GetArrayOfStructs();
    ParticleType* pstruct = aos().dataPtr();

    auto& soa = particles.GetStructOfArrays();
    auto p_realarray = soa.realarray();
    auto p_intarray = soa.intarray();

    const int phase = m_ic_solid.phase;
    const int local_cg_dem=m_is_cg_dem;

    RealBox const* const p_regions = m_p_regions;
    int const regions_nb = m_regions_nb;
    int const allow_overlap = m_allow_overlap;

    amrex::ParallelForRNG(particles_count, [pstruct,p_realarray,p_intarray,
        ic_u_s,ic_v_s,ic_w_s,mpv,phase,a_id,a_cpu,local_cg_dem,
        diameter, density,current_size,
        a_generator,p_regions,regions_nb,allow_overlap,has_granular_temperature]
      AMREX_GPU_DEVICE (int p, RandomEngine const& engine) noexcept
    {
      const int p_tot = current_size + p;

      ParticleType& part = pstruct[p_tot];

      RealVect position = a_generator.template get_position<run_on>(p, engine);

      part.pos(0) = position[0];
      part.pos(1) = position[1];
      part.pos(2) = position[2];

      part.id() = a_id+p;
      part.cpu() = a_cpu;

      if (!allow_overlap) {
        for (int box(0); box < regions_nb-1; box++) {
          if (p_regions[box].contains(part.pos())) {
            part.id().make_invalid();
          }
        }
      }

      if (part.id() >= 0) {

        Real dp = diameter.template sample<run_on>(p, engine);

        // Compute the weight using the "multi particle volume" and
        // the "single particle" diameter.
        Real const statwt = ( mpv == 0. ) ? 1.0 : mpv / ( (M_PI/6.0)*dp*dp*dp);

        // For cg-dem, scale the particle diameter using the stastical weight.
        // This is diameter is what is used to compute collisions.
        if (local_cg_dem) { dp *= std::cbrt(statwt); }

        Real rad = 0.5*dp;

        Real rho = density.template sample<run_on>(engine);

        if (has_granular_temperature) {
          p_realarray[SoArealData::velx][p_tot] = amrex::RandomNormal(0., 1., engine);
          p_realarray[SoArealData::vely][p_tot] = amrex::RandomNormal(0., 1., engine);
          p_realarray[SoArealData::velz][p_tot] = amrex::RandomNormal(0., 1., engine);
        } else {
          p_realarray[SoArealData::velx][p_tot] = ic_u_s;
          p_realarray[SoArealData::vely][p_tot] = ic_v_s;
          p_realarray[SoArealData::velz][p_tot] = ic_w_s;
        }


        p_realarray[SoArealData::statwt][p_tot] = statwt;

        p_realarray[SoArealData::radius][p_tot] = rad;
        p_realarray[SoArealData::density][p_tot] = rho;

        p_realarray[SoArealData::omegax][p_tot] = 0.0;
        p_realarray[SoArealData::omegay][p_tot] = 0.0;
        p_realarray[SoArealData::omegaz][p_tot] = 0.0;

        p_realarray[SoArealData::drag_coeff][p_tot] = 0.0;

        p_realarray[SoArealData::vel_source_x][p_tot] = 0.0;
        p_realarray[SoArealData::vel_source_y][p_tot] = 0.0;
        p_realarray[SoArealData::vel_source_z][p_tot] = 0.0;

        p_intarray[SoAintData::phase][p_tot] = phase;
        p_intarray[SoAintData::state][p_tot] = state::normal;
#if MFIX_POLYDISPERSE
        p_intarray[SoAintData::ptype][p_tot] = 0;
#endif
      }
    });
  }

  return particles_count;
}
