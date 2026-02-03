#include <restarter.H>
#include <deposition/mfix_deposition_K.H>
#include <mfix_indexes_aux.H>

#include <AMReX_ParmParse.H>

using namespace amrex;

MFIXRestarter::MFIXRestarter (const int nlev_in)
  : m_refinement_ratio(1)
  , m_eps_tolerance(1.e-15)
  , m_eps_overflow(1.)
  , m_inputs_pdiameter(0.)
  , m_inputs_pdensity(0.)
  , m_add_thermal_noise(0)
  , nlev(nlev_in)
  , avgdPIC_coarse(nlev_in, nullptr)
  , avgdPIC_fine(nlev_in, nullptr)
{
  ParmParse pp("pic2dem");

  pp.query("refinement_ratio", m_refinement_ratio);
  pp.query("eps_tolerance", m_eps_tolerance);
  pp.query("eps_overflow", m_eps_overflow);
  pp.query("thermal_noise", m_add_thermal_noise);

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_refinement_ratio > 0,
      "Error: refinement ratio must be a positive integer");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_eps_tolerance > 0.,
      "Error: eps tolerance must be a positive real number");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_eps_overflow > 0,
      "Error: eps overflow must be a positive real number");
}


MFIXRestarter::~MFIXRestarter ()
{
  for (int lev(0); lev < nlev; ++lev) {
    if (avgdPIC_coarse[lev] != nullptr)
      delete avgdPIC_coarse[lev];

    if (avgdPIC_fine[lev] != nullptr)
      delete avgdPIC_fine[lev];
  }
}


void
MFIXRestarter::allocate_avgdPIC_coarse (const mfix* mfix_coarse)
{
  auto& solids = mfix_coarse->solids;
  Transfer txfr_idxs(solids);

  for (int lev(0); lev < nlev; ++lev) {
    const auto* epf_coarse = mfix_coarse->leveldata_const().epf_const(lev);

    avgdPIC_coarse[lev] = new MultiFab(epf_coarse->boxArray(), epf_coarse->DistributionMap(),
        txfr_idxs.count, epf_coarse->nGrow(), MFInfo(), epf_coarse->Factory());
  }
}


void
MFIXRestarter::allocate_avgdPIC_fine (const mfix* mfix_fine)
{
  auto& solids = mfix_fine->solids;
  Transfer txfr_idxs(solids);

  for (int lev(0); lev < nlev; ++lev) {
    const auto* epf_fine = mfix_fine->leveldata_const().epf_const(lev);

    avgdPIC_fine[lev] = new MultiFab(epf_fine->boxArray(), epf_fine->DistributionMap(),
        txfr_idxs.count, epf_fine->nGrow(), MFInfo(), epf_fine->Factory());
  }
}


void
MFIXRestarter::change_inputs_table ()
{
  ParmParse pp_eb2("eb2");
  ParmParse pp_mfix("mfix");
  ParmParse pp_pic2dem("pic2dem");

  // Small volfrac
  {
    Real small_volfrac(0.);
    int contains_small_volfrac = pp_pic2dem.query("small_volfrac", small_volfrac);

    if (contains_small_volfrac) {
      if (pp_eb2.contains("small_volfrac"))
        pp_eb2.remove("small_volfrac");

      pp_eb2.add("small_volfrac", small_volfrac);
    }
  }

  // Geom chk file
  {
    if (pp_mfix.contains("geom_chk_file"))
      pp_mfix.remove("geom_chk_file");

    std::string geom_chk_file("");
    if (pp_pic2dem.query("geom_chk_file", geom_chk_file))
      pp_mfix.add("geom_chk_file", geom_chk_file);
  }

  // Geom levelset chk file
  {
    if (pp_mfix.contains("geom_levelset_chk_file"))
      pp_mfix.remove("geom_levelset_chk_file");

    std::string geom_levelset_chk_file("");
    if (pp_pic2dem.query("geom_levelset_chk_file", geom_levelset_chk_file))
      pp_mfix.add("geom_levelset_chk_file", geom_levelset_chk_file);
  }

  // Geom chk write
  {
    if (pp_mfix.contains("geom_chk_write"))
      pp_mfix.remove("geom_chk_write");

    bool geom_chk_write(0);
    if (pp_pic2dem.query("geom_chk_write", geom_chk_write))
      pp_mfix.add("geom_chk_write", geom_chk_write);
  }

  // Geom chk read
  {
    if (pp_mfix.contains("geom_chk_read"))
      pp_mfix.remove("geom_chk_read");

    bool geom_chk_read(0);
    if (pp_pic2dem.query("geom_chk_read", geom_chk_read))
      pp_mfix.add("geom_chk_read", geom_chk_read);
  }

  // Geometry filename
  if (!pp_mfix.contains("geom_chk_read")) {

    std::string geometry_filename("");

    if (pp_pic2dem.query("geometry_filename", geometry_filename) &&
        !pp_mfix.contains("geometry_filename")) {

      pp_pic2dem.remove("geometry_filename");

      pp_mfix.add("geometry_filename", geometry_filename);

      if (pp_mfix.contains("geom_chk_read"))
        pp_mfix.remove("geom_chk_read");
    }

  } else if (pp_mfix.contains("geometry_filename")) {

    pp_mfix.remove("geometry_filename");

  }

  {
    // Need to set
    // - do_initial_proj = true
    // - m_initial_iterations = 0
    if (pp_mfix.contains("do_initial_proj")) pp_mfix.remove("do_initial_proj");
    if (pp_mfix.contains("initial_iterations")) pp_mfix.remove("initial_iterations");

    pp_mfix.add("do_initial_proj", true);
    pp_mfix.add("initial_iterations", 0);
  }
}


void
MFIXRestarter::
set_fine_grids_from_coarse (mfix* mfix_fine,
                            const mfix* mfix_coarse) const
{
  Vector<Geometry> geom(nlev, Geometry());
  Vector<BoxArray> ba(nlev, BoxArray());

  for (int lev(0); lev < nlev; ++lev) {
    geom[lev] = mfix_coarse->Geom(lev);
    ba[lev] = mfix_coarse->boxArray(lev);

    geom[lev].refine(IntVect(m_refinement_ratio));
    ba[lev].refine(m_refinement_ratio);

    mfix_fine->SetGeometry(lev, geom[lev]);
    mfix_fine->SetBoxArray(lev, ba[lev]);
    mfix_fine->SetDistributionMap(lev, mfix_coarse->DistributionMap(lev));
  }
}


void
MFIXRestarter::
convert_coarse_data (const mfix* mfix_coarse,
                     mfix* mfix_fine) const
{
  Transfer txfr_idxs(mfix_coarse->solids);
  const int txfr_count = txfr_idxs.count;

  const int solve_enthalpy = mfix_coarse->fluid.solve_enthalpy();
  const int solve_species  = mfix_coarse->fluid.solve_species();

  const int nspecies = mfix_coarse->fluid.nspecies();
  const int ntracer  = mfix_coarse->fluid.ntracer();

  TxfrAuxiliary aux;

  for (int lev(0); lev < nlev; ++lev) {

    const auto& leveldata_coarse = mfix_coarse->leveldata_const();
    auto& leveldata_fine = mfix_fine->leveldata();

    const int has_tracer = leveldata_coarse.m_level_data[lev]->has_tracer();
    const int has_enthalpy = leveldata_coarse.m_level_data[lev]->has_enthalpy();
    const int has_temperature = leveldata_coarse.m_level_data[lev]->has_temperature();
    const int has_species = leveldata_coarse.m_level_data[lev]->has_species();

    const auto& flags_coarse = mfix_coarse->m_eb->Factory(lev).getMultiEBCellFlagFab();
    const auto& flags_fine   = mfix_fine->m_eb->Factory(lev).getMultiEBCellFlagFab();

    const amrex::MultiFab& volfrac_coarse = mfix_coarse->m_eb->Factory(lev).getVolFrac();

    // CC values
    for (MFIter mfi(*(leveldata_fine.epf(lev)),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& box = mfi.tilebox();

      GeometryData geom_coarse_data = mfix_coarse->Geom(lev).data();
      GeometryData geom_fine_data   = mfix_fine->Geom(lev).data();

      Array4<Real> dummy_arr;
      Array4<const Real> dummy_const_arr;

      const auto& epf_coarse_arr      = leveldata_coarse.epf_const(lev,mfi);
      const auto& ro_g_coarse_arr     = leveldata_coarse.rho_const(lev,mfi);
      const auto& vel_coarse_arr      = leveldata_coarse.vel_const(lev,mfi);
      const auto& T_g_coarse_arr      = leveldata_coarse.T_const(lev,mfi);
      const auto& X_coarse_arr        = leveldata_coarse.X_const(lev,mfi);
      const auto& grad_p_coarse_arr   = leveldata_coarse.grad_p_const(lev,mfi);
      const auto& trac_coarse_arr     = leveldata_coarse.tracer_const(lev,mfi);

      const auto& flags_coarse_arr   = flags_coarse.const_array(mfi);
      const auto& volfrac_coarse_arr = volfrac_coarse.const_array(mfi);

      const auto& epf_fine_arr      = leveldata_fine.epf(lev,mfi);
      const auto& ro_g_fine_arr     = leveldata_fine.rho(lev,mfi);
      const auto& vel_fine_arr      = leveldata_fine.vel(lev,mfi);
      const auto& T_g_fine_arr      = leveldata_fine.T(lev,mfi);
      const auto& X_fine_arr        = leveldata_fine.X(lev,mfi);
      const auto& grad_p_fine_arr   = leveldata_fine.grad_p(lev,mfi);
      const auto& trac_fine_arr     = leveldata_fine.tracer(lev,mfi);

      const auto& flags_fine_arr = flags_fine.const_array(mfi);

      const auto& avgdPIC_coarse_arr = avgdPIC_coarse[lev]->const_array(mfi);
      const auto& avgdPIC_fine_arr = avgdPIC_fine[lev]->array(mfi);

      amrex::ParallelFor(box, [geom_coarse_data,geom_fine_data,epf_coarse_arr,
          ro_g_coarse_arr,trac_coarse_arr,vel_coarse_arr,grad_p_coarse_arr,
          T_g_coarse_arr,epf_fine_arr,ro_g_fine_arr,trac_fine_arr,T_g_fine_arr,
          vel_fine_arr,grad_p_fine_arr,avgdPIC_coarse_arr,avgdPIC_fine_arr,
          txfr_count,solve_enthalpy,solve_species,X_fine_arr,X_coarse_arr,
          nspecies,ntracer,flags_coarse_arr,flags_fine_arr,volfrac_coarse_arr,aux,
          has_species, has_tracer, has_enthalpy, has_temperature ]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        // Create IntVect from fine indexes
        IntVect ijk_fine(i,j,k);

        if (!flags_fine_arr(i,j,k).isCovered()) {

          // Get coordinates corresponding to fine indexes
          RealVect coords = aux.get_coordinates_cc(geom_fine_data, ijk_fine);

          // Get coarse indexes
          IntVect ijk_coarse = aux.get_indexes_cc(geom_coarse_data, coords);

          GpuArray<GpuArray<GpuArray<Real,3>,3>,3> weights;

          const Box& domain_coarse = geom_coarse_data.domain;

          aux.get_weights_cc(weights, volfrac_coarse_arr, epf_coarse_arr,
              domain_coarse, ijk_coarse, flags_coarse_arr, ijk_fine, flags_fine_arr);

          epf_fine_arr(ijk_fine) = aux.sum_up(weights, domain_coarse,
              flags_coarse_arr, epf_coarse_arr, ijk_coarse);

          ro_g_fine_arr(ijk_fine) = aux.sum_up(weights, domain_coarse,
              flags_coarse_arr, ro_g_coarse_arr, ijk_coarse);

          for (int n(0); n < AMREX_SPACEDIM; ++n) {
            vel_fine_arr(ijk_fine,n) = aux.sum_up(weights, domain_coarse,
                flags_coarse_arr, vel_coarse_arr, ijk_coarse, n);

            grad_p_fine_arr(ijk_fine,n) = aux.sum_up(weights, domain_coarse,
                flags_coarse_arr, grad_p_coarse_arr, ijk_coarse, n);
          }

          if (has_tracer) {
            for (int n(0); n < ntracer; ++n) {
              trac_fine_arr(ijk_fine) = aux.sum_up(weights, domain_coarse,
                  flags_coarse_arr, trac_coarse_arr, ijk_coarse);
            }
          }


          if (has_temperature) {
            T_g_fine_arr(ijk_fine) = aux.sum_up(weights, domain_coarse,
                flags_coarse_arr, T_g_coarse_arr, ijk_coarse);
          }

          if (has_species) {
            Real sum(0.);

            for (int n(0); n < nspecies; ++n) {
              X_fine_arr(ijk_fine,n) = aux.sum_up(weights, domain_coarse,
                  flags_coarse_arr, X_coarse_arr, ijk_coarse, n);
              sum += X_fine_arr(ijk_fine,n);
            }

            AMREX_ALWAYS_ASSERT(sum > 1.e-15);

            Real new_sum(0.);

            for (int n(0); n < nspecies; ++n) {
              X_fine_arr(ijk_fine,n) /= sum;
              new_sum += X_fine_arr(ijk_fine,n);
            }

            AMREX_ALWAYS_ASSERT(std::abs(new_sum-1.) < 1.e-15);
          }

          for (int n(0); n < txfr_count; ++n) {
            avgdPIC_fine_arr(ijk_fine,n) = aux.sum_up(weights, domain_coarse,
                flags_coarse_arr, avgdPIC_coarse_arr, ijk_coarse, n);
          }

        } else {

          epf_fine_arr(ijk_fine) = mfix::covered_val;
          ro_g_fine_arr(ijk_fine) = mfix::covered_val;

          for (int n(0); n < AMREX_SPACEDIM; ++n) {
            vel_fine_arr(ijk_fine,n) = mfix::covered_val;
            grad_p_fine_arr(ijk_fine,n) = mfix::covered_val;
          }

          if (has_tracer) {
            trac_fine_arr(ijk_fine) = mfix::covered_val;
          }

          if (has_temperature) {
            T_g_fine_arr(ijk_fine) = mfix::covered_val;
          }

          if (has_species) {
            for (int n(0); n < nspecies; ++n)
              X_fine_arr(ijk_fine,n) = mfix::covered_val;
          }

          for (int n(0); n < txfr_count; ++n) {
            avgdPIC_fine_arr(ijk_fine,n) = mfix::covered_val;
          }
        }
      });
    }

    // Nodal values
    for (MFIter mfi(*(leveldata_fine.pert_p(lev)),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box box = mfi.tilebox();

      GeometryData geom_coarse_data = mfix_coarse->Geom(lev).data();
      GeometryData geom_fine_data   = mfix_fine->Geom(lev).data();

      Array4<Real> dummy_arr;
      Array4<const Real> dummy_const_arr;

      const auto& p_g_coarse_arr = leveldata_coarse.pert_p(lev,mfi);
      const auto& p_g_fine_arr = leveldata_fine.pert_p(lev,mfi);

      amrex::ParallelFor(box, [geom_coarse_data,geom_fine_data,p_g_coarse_arr,
          p_g_fine_arr,aux]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        // Create IntVect from fine indexes
        IntVect ijk_fine(i,j,k);

        // Get coordinates corresponding to fine indexes
        RealVect coords = aux.get_coordinates_nd(geom_fine_data, ijk_fine);

        // Get coarse indexes
        IntVect ijk_coarse = aux.get_indexes_nd(geom_coarse_data, coords);

        p_g_fine_arr(ijk_fine) = p_g_coarse_arr(ijk_coarse);

        if (p_g_fine_arr(ijk_fine) > 1.e35) {
          for (int ii(-1); ii <= 1; ++ii)
          for (int jj(-1); jj <= 1; ++jj)
          for (int kk(-1); kk <= 1; ++kk) {
            IntVect ijk_shift = ijk_coarse+IntVect(ii,jj,kk);
            if (geom_coarse_data.Domain().contains(ijk_shift)) {
              printf("p_g_coarse(%d,%d,%d) = %e\n", ijk_shift[0], ijk_shift[1],
                  ijk_shift[2], p_g_coarse_arr(ijk_shift));
            }
          }
        }
      });
    }
  }
}


void
MFIXRestarter::deposit_PIC (const mfix* mfix_coarse)
{
  for (int lev = 0; lev < nlev; lev++) {
    avgdPIC_coarse[lev]->setVal(0);
  }

  {
    BCList const& bc_list = mfix_coarse->get_bc_list();

    MFIXDepOpPIC2DEM pic2demDepOp(bc_list, mfix_coarse->pc,
      &(mfix_coarse->m_eb->ParticleFactory(0)), avgdPIC_coarse);


    // Deposit the interphase transfer forces to the grid and reduce to level-0.
    pic2demDepOp.Deposit();
    pic2demDepOp.Redistribute();
    pic2demDepOp.CopyToDest();
  }

  // Divide
  for (int lev = 0; lev < nlev; lev++) {

    const auto& solids = mfix_coarse->solids;

    Transfer txfr_idxs(solids);
    const int idx_eps = txfr_idxs.idx_eps;
    const int idx_density = txfr_idxs.idx_density;
    const int idx_vel = txfr_idxs.idx_vel;
    const int idx_temp = txfr_idxs.idx_temp;
    const int idx_species = txfr_idxs.idx_species;

    const int nspecies_s = solids.nspecies();

    EB_set_covered(*avgdPIC_coarse[lev], mfix::covered_val);

    for (MFIter mfi(*avgdPIC_coarse[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& box = mfi.tilebox();

      const auto& avgdPIC_coarse_arr = avgdPIC_coarse[lev]->array(mfi);

      const int solve_enthalpy = solids.solve_enthalpy();
      const int solve_species  = solids.solve_species();

      amrex::ParallelFor(box, [avgdPIC_coarse_arr,solve_enthalpy,solve_species,idx_eps,
          idx_vel,idx_species,idx_temp,nspecies_s,idx_density]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real eps = avgdPIC_coarse_arr(i,j,k,idx_eps);

        if (eps < 1.e-15) {

          avgdPIC_coarse_arr(i,j,k,idx_density) = 0.;

          avgdPIC_coarse_arr(i,j,k,idx_vel+0) = 0.;
          avgdPIC_coarse_arr(i,j,k,idx_vel+1) = 0.;
          avgdPIC_coarse_arr(i,j,k,idx_vel+2) = 0.;

          if (solve_enthalpy)
            avgdPIC_coarse_arr(i,j,k,idx_temp) = 0.;

          if (solve_species)
            for (int n_s(0); n_s < nspecies_s; ++n_s)
              avgdPIC_coarse_arr(i,j,k,idx_species+n_s) = 0.;

        } else {

          avgdPIC_coarse_arr(i,j,k,idx_density) /= eps;

          avgdPIC_coarse_arr(i,j,k,idx_vel+0) /= eps;
          avgdPIC_coarse_arr(i,j,k,idx_vel+1) /= eps;
          avgdPIC_coarse_arr(i,j,k,idx_vel+2) /= eps;

          if (solve_enthalpy)
            avgdPIC_coarse_arr(i,j,k,idx_temp) /= eps;

          if (solve_species) {
            Real sum(0.);

            for (int n_s(0); n_s < nspecies_s; ++n_s) {
              avgdPIC_coarse_arr(i,j,k,idx_species+n_s) /= eps;
              sum += avgdPIC_coarse_arr(i,j,k,idx_species+n_s);
            }

            AMREX_ALWAYS_ASSERT(sum > 1.e-15);

            Real new_sum = 0.;

            for (int n_s(0); n_s < nspecies_s; ++n_s) {
              avgdPIC_coarse_arr(i,j,k,idx_species+n_s) /= sum;
              new_sum += avgdPIC_coarse_arr(i,j,k,idx_species+n_s);
            }

            AMREX_ALWAYS_ASSERT(std::abs(new_sum-1.) < 1.e-15);
          }
        }
      });
    }

    EB_set_covered(*avgdPIC_coarse[lev], mfix::covered_val);
  }
}



void
MFIXRestarter::get_particles_radius (const mfix* mfix_ptr)
{
  int pdiameter_assigned(0);

  const auto& ics = mfix_ptr->m_initial_conditions.ic();

  for (int icv(0); icv < ics.size(); ++icv) {

    const auto& solids_ics = ics[icv].solids;

    for (int lcs(0); lcs < solids_ics.size(); ++lcs) {

      const auto& solids = solids_ics[lcs];

      if (solids.volfrac > 1.e-15) {

        if (!pdiameter_assigned) {
          m_inputs_pdiameter = solids.diameter.get_mean();
          pdiameter_assigned = 1;
        } else {
          AMREX_ALWAYS_ASSERT(std::abs(m_inputs_pdiameter - solids.diameter.get_mean()) < 1.e-15);
        }
      }
    }
  }
}


void
MFIXRestarter::get_eps_coarse (const mfix* mfix_coarse)
{
  int lev = 0;

  const auto& solids = mfix_coarse->solids;

  Transfer txfr_idxs(solids);
  const int idx_eps = txfr_idxs.idx_eps;

  const auto& geom_coarse = mfix_coarse->Geom(lev);

  const auto& f_ba = mfix_coarse->boxArray(lev);
  const auto& f_dmap = mfix_coarse->DistributionMap(lev);

  const BoxArray& p_ba = mfix_coarse->m_eb->ParticleFactory(lev).boxArray();
  const DistributionMapping& p_dmap = mfix_coarse->m_eb->ParticleFactory(lev).DistributionMap();

  const EBFArrayBoxFactory& f_ebfactory = mfix_coarse->m_eb->Factory(lev);
  const EBFArrayBoxFactory& p_ebfactory = mfix_coarse->m_eb->ParticleFactory(lev);

  bool OnSameGrids = ((f_dmap == p_dmap) && (f_ba.CellEqual(p_ba)));

  if (OnSameGrids) {

    avgd_eps_coarse = new MultiFab(f_ba, f_dmap, 1, 0, MFInfo(), f_ebfactory);
    avgd_eps_coarse->setVal(0.);
    MultiFab::Copy(*avgd_eps_coarse, *avgdPIC_coarse[lev], idx_eps, 0, 1, 0);

  } else {

    avgd_eps_coarse = new MultiFab(p_ba, p_dmap, 1, 0, MFInfo(), p_ebfactory);
    avgd_eps_coarse->setVal(0.);
    avgd_eps_coarse->ParallelCopy(*avgdPIC_coarse[lev], idx_eps, 0, 1, 0, 0);
  }

  avgd_eps_coarse->FillBoundary(geom_coarse.periodicity());
  EB_set_covered(*avgd_eps_coarse, mfix::covered_val);
}


void
MFIXRestarter::generate_particles (const Geometry& geom_coarse,
                                   mfix* mfix_fine) const
{
  // Allocate the particle data
  if (mfix_fine->m_dem.solve())
  {
    MFIXParticleContainer* pc_fine = mfix_fine->pc;

    Real strt_init_part = ParallelDescriptor::second();

    pc_fine->AllocData();

    amrex::Print() << "Auto generating particles ..." << std::endl;

    {
      int lev = 0;

      const auto& solids = mfix_fine->solids;

      Transfer txfr_idxs(solids);

      const RealVect dx(geom_coarse.CellSize());
      const RealVect plo(geom_coarse.ProbLo());

      Long total_np = 0;

      const Real pdiameter = m_inputs_pdiameter;

      for (MFIter mfi = pc_fine->MakeMFIter(lev,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Now that we know pcount, go ahead and create a particle container for this
        // grid and add the particles to it
        auto& particles = pc_fine->DefineAndReturnParticleTile(lev, mfi);

        // Get the fine box and coarsen it
        Box bx = mfi.tilebox();
        bx.coarsen(m_refinement_ratio);

        Gpu::DeviceVector<Hex_ClosePack> hcp_vector;
        Gpu::DeviceVector<Long> gen_indexes;
        Gpu::DeviceVector<Long> gen_number;

        const int bx_numPts = bx.numPts();

        hcp_vector.clear();
        gen_indexes.clear();
        gen_number.clear();
        hcp_vector.resize(bx_numPts, Hex_ClosePack());
        gen_indexes.resize(bx_numPts, 0);
        gen_number.resize(bx_numPts, 0);

        Hex_ClosePack* hcp_vector_ptr = hcp_vector.dataPtr();
        Long* gen_number_ptr = gen_number.dataPtr();
        Long* gen_indexes_ptr = gen_indexes.dataPtr();

        Array4<const Real> const& eps_arr = avgd_eps_coarse->const_array(mfi);

        const auto& eps_fab = static_cast<EBFArrayBox const&>((*avgd_eps_coarse)[mfi]);
        const auto& flags = eps_fab.getEBCellFlagFab();
        Array4<EBCellFlag const> const& flags_arr = flags.const_array();

        const amrex::Real eps_tolerance = m_eps_tolerance;
        const amrex::Real eps_overflow = m_eps_overflow;

        // compute nb of particles for each cell
        amrex::ParallelFor(bx, [plo,dx,bx,eps_arr,flags_arr,hcp_vector_ptr,
            gen_number_ptr,eps_overflow,eps_tolerance,pdiameter]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          amrex::Real eps = eps_arr(i,j,k)*eps_overflow;

          amrex::RealVect dlo(plo[0] + dx[0]*i, plo[1] + dx[1]*j, plo[2] + dx[2]*k);
          amrex::RealVect dhi(plo[0] + dx[0]*(i+1), plo[1] + dx[1]*(j+1), plo[2] + dx[2]*(k+1));

          const Long n = get_global_index(IntVect(i,j,k), bx);

          hcp_vector_ptr[n].initialize(plo, dx);

          RealBox region(dlo.begin(), dhi.begin());

          if (flags_arr(i,j,k).isCovered())
            eps = 1.e+40;

          amrex::Box const ijk_bx(IntVect(i,j,k), IntVect(i,j,k));
          hcp_vector_ptr[n].setup(ijk_bx, region, pdiameter, eps, eps_tolerance);

          if (eps > eps_tolerance && eps < 1.) {

            gen_number_ptr[n] = hcp_vector_ptr[n].get_particles_number();

          } else {

            AMREX_ALWAYS_ASSERT(hcp_vector_ptr[n].get_particles_number() == 0);
            gen_number_ptr[n] = 0;
          }
        });

        const int vec_size = gen_indexes.size();

        AMREX_ASSERT(vec_size == bx.numPts());

        Gpu::inclusive_scan(gen_number.begin(), gen_number.end(), gen_indexes.begin());

        Long local_np(0);

#ifdef AMREX_USE_GPU
        Gpu::HostVector<Long> gen_indexes_host(vec_size, 0);
        Gpu::copy(Gpu::deviceToHost, gen_indexes.begin(), gen_indexes.end(), gen_indexes_host.begin());

        local_np = gen_indexes_host[vec_size-1];
#else
        local_np = gen_indexes[vec_size-1];
#endif

        if (local_np > 0) {

          particles.resize(local_np);

          const int id = MFIXParticleContainer::ParticleType::NextID();
          const int cpu = ParallelDescriptor::MyProc();

          // generate positions
          generate_particles(local_np, bx, particles, hcp_vector_ptr, gen_number_ptr,
                             gen_indexes_ptr, id, cpu);

          // Update the particles NextID
          MFIXParticleContainer::ParticleType::NextID(id+local_np);

          // Add components for each of the runtime variables
          const int start = SoArealData::count;
          for (int comp(0); comp < pc_fine->m_runtimeRealData.count; ++comp)
            particles.push_back_real(start+comp, local_np, 0.);

          total_np += local_np;
        }
      }

      ParallelDescriptor::ReduceLongSum(total_np);

      amrex::Print() << "Total number of generated particles: " << total_np << std::endl;

      // We shouldn't need this if the particles are tiled with one tile per grid, but otherwise
      // we do need this to move particles from tile 0 to the correct tile.
      pc_fine->Redistribute();
    }

    pc_fine->Redistribute();

    Real end_init_part = ParallelDescriptor::second() - strt_init_part;
    ParallelDescriptor::ReduceRealMax(end_init_part, ParallelDescriptor::IOProcessorNumber());
    Print() << "Time spent in initializing particles " << end_init_part << std::endl;
  }
}


void
MFIXRestarter::
generate_particles (const Long /*particles_count*/,
                    const Box& bx,
                    MFIXParticleContainer::ParticleTileType& particles,
                    const Hex_ClosePack* hcp_vector_ptr,
                    const Long* gen_number_ptr,
                    const Long* gen_indexes_ptr,
                    const int id,
                    const int cpu) const
{
  auto ptile_data = particles.getParticleTileData();

  auto& aos = particles.GetArrayOfStructs();
  MFIXParticleContainer::ParticleType* pstruct = aos().dataPtr();

  auto& soa = particles.GetStructOfArrays();
  auto p_realarray = soa.realarray();
  auto p_intarray = soa.intarray();

  amrex::ParallelForRNG(bx, [pstruct,p_realarray,p_intarray,id,cpu,ptile_data,
      hcp_vector_ptr,gen_number_ptr,gen_indexes_ptr,bx]
    AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
  {
    const Long n = get_global_index(IntVect(i,j,k), bx);

    const int np = gen_number_ptr[n];

    const Long glob_id = gen_indexes_ptr[n] - np;

    for (int nn = 0; nn < np; ++nn) {

      const Long p_tot = glob_id + nn;

      MFIXParticleContainer::ParticleType& part = pstruct[p_tot];

      part.id() = id + p_tot;
      part.cpu() = cpu;

      RealVect position = hcp_vector_ptr[n].template get_position<run_on>(nn, engine);

      part.pos(0) = position[0];
      part.pos(1) = position[1];
      part.pos(2) = position[2];

      p_realarray[SoArealData::velx][p_tot] = 9.87654321e32;
      p_realarray[SoArealData::vely][p_tot] = 9.87654321e32;
      p_realarray[SoArealData::velz][p_tot] = 9.87654321e32;

      p_realarray[SoArealData::statwt][p_tot] = 1;

      p_realarray[SoArealData::radius][p_tot] = 0.0;
      p_realarray[SoArealData::density][p_tot] = 0.0;

      p_realarray[SoArealData::omegax][p_tot] = 0.0;
      p_realarray[SoArealData::omegay][p_tot] = 0.0;
      p_realarray[SoArealData::omegaz][p_tot] = 0.0;

      p_realarray[SoArealData::drag_coeff][p_tot] = 0.0;

      p_realarray[SoArealData::vel_source_x][p_tot] = 0.0;
      p_realarray[SoArealData::vel_source_y][p_tot] = 0.0;
      p_realarray[SoArealData::vel_source_z][p_tot] = 0.0;

      p_intarray[SoAintData::phase][p_tot] = 1;
      p_intarray[SoAintData::state][p_tot] = 1;
    }
  });

  return;
}


void
MFIXRestarter::
init_particles_data (mfix* mfix_fine) const
{
  MFIXParticleContainer* pc_fine = mfix_fine->pc;
  const auto& solids = mfix_fine->solids;
  const auto& fluid_parms = mfix_fine->fluid.parameters<run_on>();
  const auto fluid_props = mfix_fine->fluid.props.data<run_on>();

  const int solve_species = solids.solve_species();
  const int solve_enthalpy = solids.solve_enthalpy();

  const int nspecies_s = solids.nspecies();
  const int idx_species = pc_fine->m_runtimeRealData.X_sn;

  const Real pdiameter = m_inputs_pdiameter;

  Transfer fld_transfer(solids);

  for (int lev = 0; lev < nlev; lev++) {

    auto& leveldata_fine = mfix_fine->leveldata();

    const auto& geom = mfix_fine->Geom(lev);

    const auto& f_ba = mfix_fine->boxArray(lev);
    const auto& f_dmap = mfix_fine->DistributionMap(lev);

    const BoxArray& p_ba = mfix_fine->m_eb->ParticleFactory(lev).boxArray();
    const DistributionMapping& p_dmap = mfix_fine->m_eb->ParticleFactory(lev).DistributionMap();

    MultiFab* avgdPIC_ptr(nullptr);
    MultiFab* fine_vel_ptr(nullptr);
    MultiFab* fine_ro_g_ptr(nullptr);

    bool OnSameGrids = ((f_dmap == p_dmap) && (f_ba.CellEqual(p_ba)));

    if (OnSameGrids) {

      avgdPIC_ptr = avgdPIC_fine[lev];
      fine_vel_ptr = leveldata_fine.vel(lev);
      fine_ro_g_ptr = leveldata_fine.rho(lev);

    } else {

      avgdPIC_ptr = new MultiFab(p_ba, p_dmap, avgdPIC_fine[lev]->nComp(), 0);
      avgdPIC_ptr->setVal(0.);

      avgdPIC_ptr->ParallelCopy(*avgdPIC_fine[lev], 0, 0, avgdPIC_fine[lev]->nComp(), 0, 0);

      fine_vel_ptr = new MultiFab(p_ba, p_dmap, 3, 0);
      fine_vel_ptr->setVal(0.);
      fine_vel_ptr->ParallelCopy(*leveldata_fine.vel(lev), 0, 0, 3, 0, 0);

      fine_ro_g_ptr = new MultiFab(p_ba, p_dmap, 1, 0);
      fine_ro_g_ptr->setVal(0.);
      fine_ro_g_ptr->ParallelCopy(*(leveldata_fine.rho(lev)), 0, 0, 1, 0, 0);
    }

    avgdPIC_ptr->FillBoundary(geom.periodicity());
    fine_vel_ptr->FillBoundary(geom.periodicity());
    fine_ro_g_ptr->FillBoundary(geom.periodicity());

    const auto dxi = geom.InvCellSizeArray();
    const auto plo = geom.ProbLoArray();

    for (MFIXParticleContainer::MFIXParIter pti(*pc_fine,lev); pti.isValid(); ++pti) {

      MFIXParticleContainer::PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& ptile = pc_fine->GetParticles(lev)[index];
      auto ptile_data = ptile.getParticleTileData();

      auto& particles = pti.GetArrayOfStructs();
      MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const int add_thermal_noise = m_add_thermal_noise;

      const int np = particles.size();

      const auto& avgdPIC_array = avgdPIC_ptr->const_array(pti);
      const auto& avgdPIC_X_array = avgdPIC_ptr->const_array(pti,fld_transfer.idx_species);
      const auto& vel_array = fine_vel_ptr->const_array(pti);
      const auto& ro_g_array = fine_ro_g_ptr->const_array(pti);

      ParallelForRNG(np, [pstruct,p_realarray,ptile_data,solve_species,plo,
          avgdPIC_array,avgdPIC_X_array,vel_array,ro_g_array,nspecies_s,idx_species,dxi,
          solve_enthalpy,pdiameter,add_thermal_noise,fld_transfer,fluid_parms,fluid_props]
        AMREX_GPU_DEVICE (int p, RandomEngine const& engine) noexcept
      {
        MFIXParticleContainer::ParticleType& part = pstruct[p];

        const RealVect& pos = part.pos();

        const RealVect lx((pos[0] - plo[0])*dxi[0],
                          (pos[1] - plo[1])*dxi[1],
                          (pos[2] - plo[2])*dxi[2]);

        const IntVect ijk = lx.floor();

        int i = ijk[0]; int j = ijk[1]; int k = ijk[2];

        p_realarray[SoArealData::radius][p] = .5*pdiameter;
        p_realarray[SoArealData::density][p] = avgdPIC_array(i,j,k,fld_transfer.idx_density);

        // Particle velocity
        RealVect avg_vel(avgdPIC_array(i,j,k,fld_transfer.idx_vel+0),
                         avgdPIC_array(i,j,k,fld_transfer.idx_vel+1),
                         avgdPIC_array(i,j,k,fld_transfer.idx_vel+2));

        if (add_thermal_noise) {
          // Fluid velocity
          const RealVect vel(vel_array(i,j,k,0),
                             vel_array(i,j,k,1),
                             vel_array(i,j,k,2));

          // Relative velocity
          const Real vrel = (avg_vel-vel).vectorLength();

          // Fluid and particle densities
          const Real ro_g = ro_g_array(i,j,k);
          const Real ro_p = p_realarray[SoArealData::density][p];

          // Fluid viscosity
          const Real T_p = solve_enthalpy ? avgdPIC_array(i,j,k,fld_transfer.idx_temp) : 0;
          const Real mu_g = fluid_props.molViscosity(T_p, IntVect(i,j,k), avgdPIC_X_array);

          // Reynolds number
          const Real Re = (pdiameter * vrel * ro_g) / mu_g;

          // Square root of granular temperature
          const Real sqrt_theta = (mu_g * 2.108 * std::pow(Re, 0.85)) /
            (pdiameter * std::sqrt(ro_g*ro_p));

          avg_vel[0] += sqrt_theta * amrex::RandomNormal(0, 1, engine);
          avg_vel[1] += sqrt_theta * amrex::RandomNormal(0, 1, engine);
          avg_vel[2] += sqrt_theta * amrex::RandomNormal(0, 1, engine);
        }

        p_realarray[SoArealData::velx][p] = avg_vel[0];
        p_realarray[SoArealData::vely][p] = avg_vel[1];
        p_realarray[SoArealData::velz][p] = avg_vel[2];

        // Add temperature and species initialization
        if (solve_enthalpy)
          p_realarray[SoArealData::temperature][p] = avgdPIC_array(i,j,k,fld_transfer.idx_temp);

        if (solve_species) {
          Real X_sum(0);

          for (int n_s(0); n_s < nspecies_s; ++n_s) {
            ptile_data.m_runtime_rdata[idx_species+n_s][p] = avgdPIC_array(i,j,k,fld_transfer.idx_species+n_s);
            X_sum += ptile_data.m_runtime_rdata[idx_species+n_s][p];
          }

          AMREX_ALWAYS_ASSERT(X_sum > 1.e-15);

          for (int n_s(0); n_s < nspecies_s; ++n_s) {
            ptile_data.m_runtime_rdata[idx_species+n_s][p] /= X_sum;
          }
        }
      });
    } // ParIter loop

    if (!OnSameGrids) {
      delete avgdPIC_ptr;
      delete fine_vel_ptr;
      delete fine_ro_g_ptr;
    }

  } // lev loop

  return;
}
