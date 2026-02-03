#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_AmrCore.H>
#include <AMReX_buildInfo.H>
#include <AMReX_Geometry.H>

#include <mfix.H>
#include <mfix_fluid.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_rw.H>

#include <mfix_versions.H>
#include <mfix_restart.H>


void
mfix::Restart (std::string& restart_file,
               int& nstep,
               Real& dt,
               Real& time,
               IntVect& Nrep)
{
  BL_PROFILE("mfix::Restart()");

  Print() << "  Restarting from checkpoint " << restart_file << '\n';

  std::string const level_prefix("Level_");

  bool const has_particles = (m_dem.solve() || m_pic.solve());
  bool const replicate = (Nrep != IntVect::TheUnitVector());

  if (ooo_debug) { Print() << "Restart\n"; }

  if (replicate) {
    Print() << "  Replication " << Nrep << '\n';

    Print() << "Due to replication, level-set will be re-calculated.\n";
    levelset_restart = false;
  }

  RealVect prob_lo;
  RealVect prob_hi;
  IntVect n_cell;

  /***************************************************************************
   * Load header: set up problem domain (including BoxArray)                 *
   *              allocate leveldata and load fluid data                     *
   *              load particle data                                         *
   ***************************************************************************/

  Real version_nb(0.);

  //{
    std::string File(restart_file + "/Header");
    std::string particlesFile(restart_file + "/particles/Header");

    int num_real_comps(0);
    int num_real_comps_read(0);

    if (has_particles) {

      // Read num_real_comps_read
      std::ifstream ff;
      ff.open(particlesFile);

      if (ff.fail()) {
         amrex::Abort("Failed to open particles checkpoint file");
      }

      std::string temp_str;
      ff >> temp_str;

      int temp_int;
      ff >> temp_int;

      ff >> num_real_comps_read;
      ff.close();
    }

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    // Auxiliary strings
    std::string line, word;

    // Version information
    std::string version_str;

    is >> version_str; AMREX_ALWAYS_ASSERT(version_str == "Checkpoint");
    is >> version_str; AMREX_ALWAYS_ASSERT(version_str == "version:");
    is >> version_str;

    version::checkpoint version(version_str);

    if ( version() == std::tuple<int,int>(1,0)) {
      reporter::Log(reporter::Warning,__FILE__, __LINE__)
          << "Can't check whether inputs match Checkpoint file or not.\n"
          << "Use checkpoint file version > 1.0";
    }

    // Multi-level control
    std::getline(is, line);
    int chkpoint_nlev;
    is >> chkpoint_nlev;
    finest_level = chkpoint_nlev-1;

    // Time stepping controls
    m_rw->GotoNextLine(is);
    is >> nstep;

    // Current time step
    m_rw->GotoNextLine(is);
    is >> dt;

    // Current time.
    m_rw->GotoNextLine(is);
    is >> time;

    // Geometry controls
    if ( version() == std::tuple<int,int>(1,0)) {

      m_rw->GotoNextLine(is);
      std::getline(is, line);
      {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
          prob_lo[i++] = std::stod(word);
        }
      }

      std::getline(is, line);
      {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
          prob_hi[i++] = std::stod(word);
        }
      }

    } else {

      m_rw->GotoNextLine(is);
      is >> prob_lo;

      m_rw->GotoNextLine(is);
      is >> prob_hi;

      m_rw->GotoNextLine(is);
      is >> n_cell;

      // Check that inputs geometry matches checkpoint geometry
      RealVect inputs_prob_lo(Geom(0).ProbLo());
      if ( prob_lo != inputs_prob_lo) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Checkpoint file prob_lo does not match inputs value!";
      }

      RealVect inputs_prob_hi(Geom(0).ProbHi());
      if ( prob_hi != inputs_prob_hi) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Checkpoint file prob_hi does not match inputs value!";
      }

      IntVect inputs_n_cell(Geom(0).Domain().size());
      if ( n_cell != inputs_n_cell) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Checkpoint file n_cel does not match inputs value!";
      }
    }

    if ( version() > std::tuple<int,int>(1,0)) {

      m_rw->GotoNextLine(is);
      Real small_volfrac(0.);
      is >> small_volfrac;

      Real inputs_small_volfrac(0.);
      ParmParse pp("eb2");
      pp.query("small_volfrac", inputs_small_volfrac);

      Real error = Math::abs(small_volfrac - inputs_small_volfrac);
      Real tolerance = std::numeric_limits<Real>::epsilon();
      if ( error > tolerance ) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Checkpoint file small_volfrac does not match inputs value!"
            << "\n  checkpoint value: " << small_volfrac
            << "\n  intpus value:     " << inputs_small_volfrac;
      }
    }


    // All the header info has been read at this point.
    // Allocate leveldata MultiFabs and copy data into them.


    // Replicate
    if (replicate) {
       for (int dir(0); dir < AMREX_SPACEDIM; dir++) {
          prob_lo[dir] *= Nrep[dir];
          prob_hi[dir] *= Nrep[dir];
          if ( version() == std::tuple<int,int>(1,1)) { n_cell[dir] *= Nrep[dir]; }
       }
    }


    // BoxArray controls
    for (int lev(0); lev < chkpoint_nlev; ++lev) {

      BoxArray orig_ba, ba;

      orig_ba.readFrom(is);
      m_rw->GotoNextLine(is);

      Box orig_domain(orig_ba.minimalBox());

      if (replicate) {
         Print() << " OLD BA had " << orig_ba.size()  << " GRIDS\n"
             << " OLD Domain" << orig_domain << '\n';
      }

      BoxList bl;
      for (int nb = 0; nb < orig_ba.size(); nb++) {
        for (int k = 0; k < Nrep[2]; k++) {
          for (int j = 0; j < Nrep[1]; j++) {
            for (int i = 0; i < Nrep[0]; i++) {
              Box b(orig_ba[nb]);
              IntVect shift_vec(i*orig_domain.length(0),
                                j*orig_domain.length(1),
                                k*orig_domain.length(2));
              b.shift(shift_vec);
              bl.push_back(b);
            }
          }
        }
      }

      ba.define(bl);
      DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

      if (replicate) {

        RealBox rb(prob_lo.begin(), prob_hi.begin());
        Geom(lev).ProbDomain(rb);
        Geom(lev).ResetDefaultProbDomain(rb);

        Box new_domain(ba.minimalBox());
        geom[lev].Domain(new_domain);
      }

      MakeNewLevelFromScratch(lev, time, ba, dm);

    } // lev

    Print() << "  Finished reading header\n";
  //}


  /***************************************************************************
   * Load fluid data                                                         *
   ***************************************************************************/
  if (fluid.solve()) {

    // Load the field data
    for (int lev(0); lev < chkpoint_nlev; ++lev) {

      Print() << "  Loading vel\n";
      restart_leveldata_var ( Nrep, leveldata().vel(lev),
          MultiFabFileFullPrefix(lev, restart_file, level_prefix, "u_g") );

      Print() << "  Loading grad_p\n";
      restart_leveldata_var ( Nrep, leveldata().grad_p(lev),
          MultiFabFileFullPrefix(lev, restart_file, level_prefix, "gpx") );

      Print() << "  Loading epf\n";
      restart_leveldata_var ( Nrep, leveldata().epf(lev),
          MultiFabFileFullPrefix(lev, restart_file, level_prefix, "ep_g") );

      Print() << "  Loading pert_p\n";
      restart_leveldata_var ( Nrep, leveldata().pert_p(lev),
          MultiFabFileFullPrefix(lev, restart_file, level_prefix, "p_g") );

      Print() << "  Loading rho_g\n";
      restart_leveldata_var ( Nrep, leveldata().rho(lev),
          MultiFabFileFullPrefix(lev, restart_file, level_prefix, "ro_g") );

      if (leveldata(lev)->has_species()) {

        Print() << "  Loading X\n";
        restart_leveldata_var ( Nrep, leveldata().X(lev),
            MultiFabFileFullPrefix(lev, restart_file, level_prefix, "X_gk") );
      }

      if (leveldata(lev)->has_temperature()) {

        Print() << "  Loading T\n";
        restart_leveldata_var ( Nrep, leveldata().T(lev),
            MultiFabFileFullPrefix(lev, restart_file, level_prefix, "T_g") );

        if (fluid.solve_enthalpy()) {

          Print() << "  Computing h from T\n";

          const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>
              (leveldata().h(lev)->Factory());

          for (MFIter mfi(*leveldata().h(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            Box const& bx = mfi.tilebox();

            EBCellFlagFab const& flagfab = factory.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flags = flagfab.const_array();

            Array4<Real const> const& X = leveldata().X_const(lev,mfi);
            Array4<Real const> const& T = leveldata().T_const(lev, mfi);

            Array4<Real      > const& h = leveldata().h(lev,mfi);

            const auto fluid_props = fluid.props.data<run_on>();

            ParallelFor(bx, [T,h,X,fluid_props,flags]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              const int cell_is_covered = static_cast<int>(flags(i,j,k).isCovered());

              h(i,j,k) = fluid_props.enthalpy(IntVect(i,j,k), T, X, cell_is_covered);

            });
          } // MFIter
        } // solve_enthalpy
      } // has temperatuer


      // Thermodynamic pressure
      if (fluid.solve_enthalpy() &&
          fluid.constraint.isIdealGasClosedSystem() ) {

        if (version_nb > 1.1) {
          // Read thermodynamic pressure from input file
          auto prefix = amrex::MultiFabFileFullPrefix(lev,
              restart_file, level_prefix, "thermodynamic_p_g");

          std::ifstream filestream;
          filestream.open(prefix, std::ios::in);

          filestream >> therm_p;

          filestream.close();

        } else {

          AMREX_ALWAYS_ASSERT(fluid.solve_species());

          amrex::Real thermo_p_g(0.);
          average_thermodynamic_pressure(lev, thermo_p_g);
          fluid.set_thermodynamic_pressure(thermo_p_g);
        }

      } else if (fluid.solve_enthalpy() && fluid.solve_species() &&
          fluid.constraint.isIdealGasOpenSystem() ) {

          amrex::Real thermo_p_g(0.);
          average_thermodynamic_pressure(lev, thermo_p_g);

          const Real rel_tol = 1.e-3;
          const Real abs_tol = 1.e-1;

          const Real IC_thermo_p_g = fluid.thermodynamic_pressure();

          if (std::abs(thermo_p_g - IC_thermo_p_g) > rel_tol*IC_thermo_p_g + abs_tol) {
            reporter::Log(reporter::Warning,__FILE__, __LINE__)
                << "Averaged thermodynamic pressure differs from inputs value!"
                << "\n  Computed average: " << thermo_p_g
                << "\n  Inputs value:     " << IC_thermo_p_g;
          }
      }

    } // lev

    Print() << "  Finished reading fluid data\n";
  } // solve fluid


  // Particle data is loaded into the MFIXParticleContainer's base
  // class using amrex::NeighborParticleContainer::Restart
  // at this step, pc is constructed from mfix.init and it may have
  // different ba and dm from the fluid grids.

  if (has_particles) {

    pc = new MFIXParticleContainer(geom, dmap[0], grids[0],
        ics(), bcs(), solids, m_dem, m_pic, fluid, reactions,
        m_coupling.include_virtual_mass());

    for ( int lev(0); lev<nlev(); ++lev) {
      m_eb->make_particle_factory( lev, geom[lev], grids[lev], dmap[lev]);
    }

    const int runtimeReal_count = pc->m_runtimeRealData.count;
    const int runtimeInt_count = pc->m_runtimeIntData.count;

    int const ncomp_ep_s    = pc->m_runtimeRealData.ncomp_eps;
    int const ncomp_energy  = pc->m_runtimeRealData.ncomp_energy;

    if ( version.is_current() ) {

      pc->Restart(restart_file, "particles");

    } else {

      num_real_comps = SoArealData::count + pc->m_runtimeRealData.count;

      Print() << "Patching particle restart"
          << "\n  Number of reals in pc: " << num_real_comps
          << "\n  Number of reals read:  " << num_real_comps_read << '\n';

      AMREX_ALWAYS_ASSERT(num_real_comps < num_real_comps_read);

      if ( version() == std::tuple<int,int>(1,5) ) {

        pcPatch_v_1_5 patch(version(), runtimeReal_count,
            runtimeInt_count, ncomp_ep_s, ncomp_energy);
        restart_from_old_chkpt(restart_file, "particles", patch);

      } else if ( version() == std::tuple<int,int>(1,4) ) {

        pcPatch_v_1_4 patch(version(), runtimeReal_count,
            runtimeInt_count, ncomp_ep_s, ncomp_energy);
        restart_from_old_chkpt(restart_file, "particles", patch);

      } else if ( version() == std::tuple<int,int>(1,3) ) {

        pcPatch_v_1_3 patch(version(), runtimeReal_count,
            runtimeInt_count, ncomp_ep_s, ncomp_energy);
        restart_from_old_chkpt(restart_file, "particles", patch);

      } else if ( version() <= std::tuple<int,int>{1,2} ) {

        pcPatch_v_1_2 patch(version(), runtimeReal_count,
            runtimeInt_count, ncomp_ep_s, ncomp_energy);
        restart_from_old_chkpt(restart_file, "particles", patch);

      } else {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Checkpoint file version is unknown!\n"
          << "It is not possible to restart the simulation.";
      }
    }

    Print() << "Finished reading particle data\n";


    // Make sure that the particle BoxArray is the same as the mesh data -- we can
    //      create a dual grid decomposition in the regrid operation
    for (int lev(0); lev <= finestLevel(); lev++) {

      pc->SetParticleBoxArray       (lev, grids[lev]);
      pc->SetParticleDistributionMap(lev,  dmap[lev]);
      pc->SetParticleGeometry       (lev,  geom[lev]);
    }

    pc->Redistribute(0, 0, 0, 0, false);

    int lev = 0;
    if (Nrep != IntVect::TheUnitVector()) {
      pc->Replicate(Nrep, geom[lev], dmap[lev], grids[lev]);
    }

    // load level set parameters from checkpoint
    if (levelset_restart) {

      int levelset_params[4] = { m_eb->levelset_refinement(),
                                 m_eb->levelset_pad(),
                                 m_eb->levelset_eb_refinement(),
                                 m_eb->levelset_eb_pad()         };

      std::ifstream param_file;
      std::stringstream param_file_name;
      param_file_name << restart_file << "/LSFactory_params";
      param_file.open(param_file_name.str());

      readIntData(levelset_params, 4, param_file, FPC::NativeIntDescriptor());

      int ls_ref = levelset_params[0], ls_pad = levelset_params[1],
          eb_ref = levelset_params[2], eb_pad = levelset_params[3];

      Print() << "     + Loaded level-set parameters:" << std::endl
              << "       ref = " << ls_ref << "    pad = " << ls_pad
              << "    eb_ref = " << eb_ref << " eb_pad = " << eb_pad
              << std::endl;

      // Inform the user if the checkpoint parameters do not match those in the
      // inputs file. The checkpoint inputs overwrite the inputs file.
      if (ls_ref != m_eb->levelset_refinement()) {
        Print() << "     * Overwrote levelset_refinement = "
             << m_eb->levelset_refinement() << " -> " << ls_ref << '\n';
      }
      if (ls_pad != m_eb->levelset_pad()) {
        Print() << "     * Overwrote levelset_pad = "
            << m_eb->levelset_pad() << " -> " << ls_pad << '\n';
      }
      if (eb_ref != m_eb->levelset_eb_refinement()) {
        Print() << "     * Overwrote levelset_eb_refinement = "
            << m_eb->levelset_eb_refinement() << " -> " << eb_ref << '\n';
      }
      if (eb_pad != m_eb->levelset_eb_pad()) {
        Print() << "     * Overwrote levelset_eb_pad = "
            << m_eb->levelset_eb_pad() << " -> " << eb_pad << '\n';
      }
    } // levelset parameters

    m_eb->fill_levelsets(geom, pc, m_boundary_conditions, m_porous_media);

  } // has particles


  Print() << "  Done with mfix::Restart\n";
}


/*
       SoArealData changes in checkpoints
    +---------------+-----+-----+-----+-----+-----+
    | version       | 1.2 | 1.3 | 1.4 | 1.5 | 1.6 |
    +---------------+-----+-----+-----+-----+-----+
    | radius        |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | volume        |  ✓  |     |     |     |     |
    | mass          |  ✓  |  ✓  |     |     |     |
    | density       |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | ep_s          |  ✓  |  ✓  |  ✓  |     |     |
    | velx          |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | vely          |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | velz          |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | omegax        |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | omegay        |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | omegaz        |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | statwt        |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | drag_coeff,   |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | vel_source_x  |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | vel_source_y  |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | vel_source_z  |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | cp_s          |  ✓  |  ✓  |  ✓  |     |     |
    | temperature   |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    | convection    |  ✓  |  ✓  |  ✓  |  ✓  |     |
    | count         |  ✓  |  ✓  |  ✓  |  ✓  |  ✓  |
    +---------------+-----+-----+-----+-----+-----+
    | change        |  0  | -1  | -1  | -2  | -1  |
    +---------------+-----+-----+-----+-----+-----+


    +---------+-----------------------------------+-----+-----+-----+-----+
    | version | runtimeRealData change log        | 1.3 | 1.4 | 1.5 | 1.6 |
    +---------+-----------------------------------+-----+-----+-----+-----+
    |  1.2    | No change.                        |  x  |  x  |  x  |  x  |
    +---------+-----------------------------------+-----+-----+-----+-----+
    |  1.3    | No change.                        |  0  |  x  |  x  |  x  |
    +---------+-----------------------------------+-----+-----+-----+-----+
    |  1.4    | No change.                        |  0  |  0  |  x  |  x  |
    +---------+-----------------------------------+-----+-----+-----+-----+
    |  1.5    | Added 'ep_s'                      | +1  | +1  |  0  |  x  |
    +---------+-----------------------------------+-----+-----+-----+-----+
    |  1.6    | Added 'energy_sources'            | +3  | +3  | +2  |  0  |
    +---------+-----------------------------------+-----+-----+-----+-----+
*/


template <typename PatchType> void mfix::
restart_from_old_chkpt ( const std::string& a_restart_file,
                         const std::string& a_name,
                         PatchType& a_patch )
{
  int const patch_soa_has_ep_s = a_patch.soa_has_ep_s;
  int const patch_rrd_has_ep_s = a_patch.rrd_has_ep_s;

  int const patch_rrd_has_energy_src = a_patch.rrd_has_energy_src;

  auto runtimeRealData = pc->m_runtimeRealData;

  int const ncomp_species = runtimeRealData.ncomp_species;
  int const scomp_species = runtimeRealData.X_sn;

  int const ncomp_eps     = runtimeRealData.ncomp_eps;
  int const scomp_eps     = runtimeRealData.ep_s;

  int const ncomp_vm      = runtimeRealData.ncomp_vm;
  int const scomp_vm      = runtimeRealData.vm_coeff;

  int const ncomp_acc     = runtimeRealData.ncomp_acc;
  int const scomp_acc     = runtimeRealData.acceleration;

  int const ncomp_energy  = runtimeRealData.ncomp_energy;
  int const scomp_energy  = runtimeRealData.energy_source;

  int const ncomp_tan_his = 3*runtimeRealData.max_contacts;
  int const scomp_tan_his = runtimeRealData.pft_neighbor_idx;

  int const NArrayReal = static_cast<int>(PatchType::SoArealData::count);
  int const NArrayInt  = static_cast<int>(PatchType::SoAintData::count);

  using PC = amrex::NeighborParticleContainer<0,0,NArrayReal,NArrayInt>;
  using PCIter = amrex::ParIter<0,0,NArrayReal,NArrayInt>;

  PC* old_pc = new PC(GetParGDB(), 1);

  // runtime variables needed to read the "old" checkpoint file
  int const runtimeReal_count = a_patch.runtimeReal_count;
  int const runtimeInt_count = a_patch.runtimeInt_count;

  // Add Real runtime variables
  for (int n(0); n < runtimeReal_count; ++n) {
    old_pc->AddRealComp(true);
  }

  // Add Int runtime variables
  for (int n(0); n < runtimeInt_count; ++n) {
    old_pc->AddIntComp(true);
   }

  old_pc->Restart(a_restart_file, a_name);

  // Check if this is needed
  old_pc->Redistribute(0, 0, 0, 0, false);

  for (int lev(0); lev < nlev(); ++lev) {
    for (PCIter pti(*old_pc,lev); pti.isValid(); ++pti) {

      const int np = pti.numParticles();

      MFIXParticleContainer::PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& old_particles = old_pc->GetParticles(lev);
      auto& old_ptile = old_particles[index];
      auto& old_aos = old_ptile.GetArrayOfStructs();
      auto old_pstruct = old_aos().dataPtr();
      auto& old_soa = old_ptile.GetStructOfArrays();
      auto old_realarray = old_soa.realarray();
      auto old_intarray = old_soa.intarray();
      auto old_tile_data = old_ptile.getParticleTileData();

      auto& particletile = pc->DefineAndReturnParticleTile(lev, pti);
      particletile.resize(np);
      auto& particles  = pc->GetParticles(lev);
      auto& ptile = particles[index];
      auto& aos = particletile.GetArrayOfStructs();
      auto pstruct = aos().dataPtr();
      auto& soa = particletile.GetStructOfArrays();
      auto realarray = soa.realarray();
      auto intarray = soa.intarray();
      auto tile_data = ptile.getParticleTileData();

      ParallelFor(np, [pstruct,old_pstruct,old_realarray,realarray,old_intarray,
          intarray,old_tile_data,tile_data,runtimeReal_count,runtimeInt_count,
          runtimeRealData, patch_soa_has_ep_s, patch_rrd_has_ep_s, patch_rrd_has_energy_src,
          ncomp_species, scomp_species, ncomp_eps, scomp_eps, ncomp_vm, scomp_vm,
          ncomp_acc, scomp_acc, ncomp_energy, scomp_energy, ncomp_tan_his, scomp_tan_his]
      AMREX_GPU_DEVICE (int i) noexcept
      {
        auto& part = pstruct[i];
        auto& old_part = old_pstruct[i];

        part.pos(0) = old_part.pos(0);
        part.pos(1) = old_part.pos(1);
        part.pos(2) = old_part.pos(2);

        part.id() = old_part.id();
        part.cpu() = old_part.cpu();

        intarray[SoAintData::phase][i] = old_intarray[SoAintData::phase][i];
        intarray[SoAintData::state][i] = old_intarray[SoAintData::state][i];
#if MFIX_POLYDISPERSE
        intarray[SoAintData::ptype][i] = old_intarray[SoAintData::ptype][i];
#endif

        realarray[SoArealData::radius][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::radius)][i];

        realarray[SoArealData::density][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::density)][i];

        realarray[SoArealData::velx][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::velx)][i];
        realarray[SoArealData::vely][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::vely)][i];
        realarray[SoArealData::velz][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::velz)][i];

        realarray[SoArealData::omegax][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::omegax)][i];
        realarray[SoArealData::omegay][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::omegay)][i];
        realarray[SoArealData::omegaz][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::omegaz)][i];

        realarray[SoArealData::statwt][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::statwt)][i];

        realarray[SoArealData::drag_coeff][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::drag_coeff)][i];

        realarray[SoArealData::vel_source_x][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::vel_source_x)][i];
        realarray[SoArealData::vel_source_y][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::vel_source_y)][i];
        realarray[SoArealData::vel_source_z][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::vel_source_z)][i];

        realarray[SoArealData::temperature][i] =
            old_realarray[static_cast<int>(PatchType::SoArealData::temperature)][i];

        int scomp_old_rrd = 0;

        // species mass fraction :: assume if in new, in old
        for (int n(0); n<ncomp_species; ++n) {
          tile_data.m_runtime_rdata[n + scomp_species][i] =
              old_tile_data.m_runtime_rdata[n + scomp_old_rrd][i];
        }
        scomp_old_rrd += ncomp_species;


        // volume fraction :: assume if in new, in old
        // Copy fluid volume fraction at particle (PIC only)
        if ( ncomp_eps ) {

          if (patch_soa_has_ep_s) {

            tile_data.m_runtime_rdata[scomp_eps][i] =
                old_realarray[static_cast<int>(PatchType::SoArealData::ep_s)][i];

          } else if (patch_rrd_has_ep_s) {

            // Assert that if it's old runtime reals,
            // then it is in the new runtime reals
            AMREX_ASSERT(ncomp_eps == 1);

            tile_data.m_runtime_rdata[scomp_eps][i] =
                old_tile_data.m_runtime_rdata[scomp_old_rrd][i];

          } else if ( ncomp_eps > 0 ) {

            // volume fraction is in new, but not old.
            tile_data.m_runtime_rdata[scomp_eps][i] = 0.;
          }
        }

        // Even if we don't need it, account for the index shift
        if (patch_rrd_has_ep_s) { scomp_old_rrd += ncomp_eps; }


        // Assume if virtual mass is in current, it is in old
        for (int n(0); n<ncomp_vm; ++n) {
          tile_data.m_runtime_rdata[n+scomp_vm][i] =
              old_tile_data.m_runtime_rdata[n + scomp_old_rrd][i];
        }
        scomp_old_rrd += ncomp_vm;


        // Assume if acceleration is in current, it is in old
        for (int n(0); n<ncomp_acc; ++n) {
          tile_data.m_runtime_rdata[n+scomp_acc][i] =
              old_tile_data.m_runtime_rdata[n + scomp_old_rrd][i];
        }
        scomp_old_rrd += ncomp_acc;


        if (patch_rrd_has_energy_src) {

          // Assert that if it's old runtime reals,
          // then must be in the new runtime reals
          AMREX_ASSERT(ncomp_energy > 0);

          for (int n(0); n<ncomp_energy; ++n) {
            tile_data.m_runtime_rdata[n+scomp_energy][i] =
                old_tile_data.m_runtime_rdata[n + scomp_old_rrd][i];
          }
          scomp_old_rrd += ncomp_energy;

        } else {

          // If it is in the new, but not the old, zero out the new values.
          for (int n(0); n<ncomp_energy; ++n) {
            tile_data.m_runtime_rdata[n+scomp_energy][i] = 0.;
          }
        }


        // Assume if acceleration is in current, it is in old
        for (int n(0); n<ncomp_tan_his; ++n) {
          tile_data.m_runtime_rdata[n+scomp_tan_his][i] =
              old_tile_data.m_runtime_rdata[n + scomp_old_rrd][i];
        }
        scomp_old_rrd += ncomp_tan_his;

        AMREX_ASSERT( scomp_old_rrd == runtimeReal_count );

        for (int n(0); n < runtimeInt_count; ++n)
        { tile_data.m_runtime_idata[n][i] = old_tile_data.m_runtime_idata[n][i]; }

      });
    }
  }

  delete old_pc;

}



void mfix::
average_thermodynamic_pressure ( const int lev,
                                 Real& thermo_p_g) const
{
  Real local_thermo_p_g(0.);
  int local_ncells(0);

  for (MFIter mfi(*(leveldata_const().epf_const(lev)), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();
    EBCellFlagFab const& flagfab = m_eb->factory()[lev]->getMultiEBCellFlagFab()[mfi];

    if (flagfab.getType(bx) != FabType::covered) {

      ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
      ReduceData<Real,int> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      auto const& rho_array = leveldata_const().rho_const(lev,mfi);
      auto const& Tf_array = leveldata_const().T_const(lev,mfi);
      auto const& Xn_array = leveldata_const().X_const(lev,mfi);
      auto const& flags_array = flagfab.const_array();

      const int nspecies_g = fluid.nspecies();
      const auto& fluid_parms = fluid.parameters<run_on>();
      const int fluid_is_a_mixture = fluid.isMixture();

      reduce_op.eval(bx, reduce_data, [rho_array,Tf_array,Xn_array,flags_array,
          fluid_parms,fluid_is_a_mixture,nspecies_g]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        if (!flags_array(i,j,k).isCovered()) {

          Real p_g(0.);

          Real rho = rho_array(i,j,k);
          Real Tf  = Tf_array(i,j,k);
          Real R = MFIXFluidPhase::R;

          if (!fluid_is_a_mixture) {

            p_g = rho * R * Tf * (Xn_array(i,j,k,0) / fluid_parms.get_MW_g());

          } else {

            for (int n_g(0); n_g < nspecies_g; ++n_g) {
              p_g += rho * R * Tf * (Xn_array(i,j,k,n_g) / fluid_parms.get_MW_gk(n_g));
            }
          }

          return {p_g, 1};
        }

        return {0., 0};
      });

      ReduceTuple host_tuple = reduce_data.value(reduce_op);
      local_thermo_p_g = amrex::get<0>(host_tuple);
      local_ncells = amrex::get<1>(host_tuple);
    }
  }

  ParallelDescriptor::ReduceRealSum(local_thermo_p_g);
  ParallelDescriptor::ReduceIntSum(local_ncells);

  thermo_p_g = local_thermo_p_g / local_ncells;
}



void mfix::
restart_leveldata_var ( IntVect a_Nrep, MultiFab* a_leveldata_MF,
                        std::string const& a_prefix )
{

  auto replicate_data = [] (MultiFab& dst, MultiFab& src) -> void
  {
    if (src.boxArray().size() > 1) {
      amrex::Abort("Replication only works if one initial grid");
    }

    const int ncomp = src.nComp();

    FArrayBox single_fab(src.boxArray()[0], ncomp, The_Pinned_Arena());
    src.copyTo(single_fab);

#ifdef AMREX_USE_GPU
    const auto nreals = single_fab.size();
    const auto nbytes = nreals*sizeof(Real);
    FArrayBox single_fab_d(single_fab.box(), single_fab.nComp());
    Gpu::htod_memcpy(single_fab_d.dataPtr(), single_fab.dataPtr(), nbytes);
#endif

    // Copy and replicate mf into velocity
    for (MFIter mfi(dst, false); mfi.isValid(); ++mfi) {

      int ib = mfi.index();

#ifdef AMREX_USE_GPU
      dst[ib].copy<RunOn::Gpu>(single_fab_d, single_fab_d.box(), 0, mfi.validbox(), 0, ncomp);
#else
      dst[ib].copy<RunOn::Host>(single_fab, single_fab.box(), 0, mfi.validbox(), 0, ncomp);
#endif
    }
  };

  MultiFab read_MF(The_Pinned_Arena());
  VisMF::Read(read_MF, a_prefix);

  if (a_Nrep == IntVect::TheUnitVector()) {

    int const ng_to_copy(0);
    int const ncomp( a_leveldata_MF->nComp());

    a_leveldata_MF->ParallelCopy(read_MF, 0, 0, ncomp, ng_to_copy, ng_to_copy);

  } else {

    replicate_data(*a_leveldata_MF, read_MF);
  }

}
