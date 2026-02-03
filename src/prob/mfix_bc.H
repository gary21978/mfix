#ifndef _MFIX_BC_H_
#define _MFIX_BC_H_

#include <AMReX_BC_TYPES.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBMultiFabUtil.H>

#include <AMReX_Vector.H>
#include <AMReX_IArrayBox.H>

#include <AMReX_Vector.H>
#include <AMReX_LO_BCTYPES.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>

#include <mfix_eb_parms.H>
#include <mfix_regions.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

#include <bits/stdc++.h>

enum struct BCFillVar { epf, rho, vel, tracer, T, h, X, none };

static std::map<BCFillVar,std::string> BCFillVarName = {
  {BCFillVar::epf,    "epf"   }, {BCFillVar::rho,    "rho"   },
  {BCFillVar::vel,    "vel"   }, {BCFillVar::tracer, "tracer"},
  {BCFillVar::T,      "T"     }, {BCFillVar::h,      "h"     },
  {BCFillVar::X,      "X"     }, {BCFillVar::none,   "none"  }
};


// Boundary conditions list auxiliary structure
struct BCList {

  public:
    // Types
    enum {undefined=0, cover=9, pinf=10, pout=11, minf=20, eb=111};

    // Constructor
    BCList (const int a_max_levels);

    // Destructor
    ~BCList();

    void MakeBCArrays (int const a_lev, int const a_nghost,
                       const amrex::Geometry& a_geom, bool const a_debug);

  public:
    // Boundary conditions types
    amrex::Vector<amrex::IArrayBox*> bc_ilo;
    amrex::Vector<amrex::IArrayBox*> bc_ihi;
    amrex::Vector<amrex::IArrayBox*> bc_jlo;
    amrex::Vector<amrex::IArrayBox*> bc_jhi;
    amrex::Vector<amrex::IArrayBox*> bc_klo;
    amrex::Vector<amrex::IArrayBox*> bc_khi;
};


// Forward declaration
struct BC_t;


class MFIXBoundaryConditions
{
  public:
    // Constructor
    MFIXBoundaryConditions ( int const a_max_level,
                            amrex::Vector<amrex::Geometry> const& a_geom,
                            amrex::Vector<amrex::IntVect> a_ref_ratio,
                            BCList& a_bc_list, MFIXEmbeddedBoundaries& a_eb);

    void Initialize ( const MFIXRegions& regions,
                      MFIXFluidPhase& fluid,
                      MFIXSolidsPhase& solids,
                      MFIXDEM& dem,
                      MFIXPIC& pic);

    void set_bc_list ( int const a_lev,
                       int const a_nghost_bc);

    void set_bc_type (int a_lev, int a_nghost_bc,
                      MFIXFluidPhase& a_fluid );

    void set_bcrec_lo (const int lev, const int dir, const int l_type);
    void set_bcrec_hi (const int lev, const int dir, const int l_type);


    void set_pressure_values ( int const a_solve_fluid );

    void set_epf_values ( int const a_solve_fluid );

    // Compute the bc_value based on the current time

    void set_density_bc_values (amrex::Real a_time);

    void set_velocity_bc_values (amrex::Real a_time);

    void set_energy_bc_values ( amrex::Real a_time,
                                const MFIXFluidPhase& fluid);

    void set_species_bc_values (amrex::Real a_time,
                                int const a_nspecies );

    void set_tracer_bc_values (amrex::Real a_time,
                               const MFIXFluidPhase& fluid);

    // Copy the bc_value into the domain boundary cells

    void set_epf_bcs ( int const a_lev, amrex::Real const a_time,
                       amrex::MultiFab* a_epf,
                       int const a_minf_type = amrex::BCType::ext_dir);

    void set_density_bcs ( int const a_lev, amrex::Real const a_time,
                           amrex::MultiFab* a_rho);

    void set_velocity_bcs ( int const a_lev, amrex::Real const a_time,
                            amrex::MultiFab* a_vel,
                            int const a_minf_type = amrex::BCType::ext_dir);

    void set_tracer_bcs ( int const a_lev, amrex::Real const a_time,
                          const MFIXFluidPhase& fluid,
                          amrex::MultiFab* a_tracer);

    void set_temperature_bcs ( int const a_lev, amrex::Real const a_time,
                               const MFIXFluidPhase& fluid,
                               amrex::MultiFab* a_T);

    void set_enthalpy_bcs (int const a_lev, amrex::Real const a_time,
                           const MFIXFluidPhase& fluid,
                           amrex::MultiFab* a_h);

    void set_species_bcs ( int const a_lev, amrex::Real const a_time,
                           int const a_nspecies, amrex::MultiFab* Xk);

    void set_vec_bcs (int const a_lev, amrex::Real const a_time,
                      amrex::MultiFab* mf_in);


    // Copy Diriclet EB values into cut-cells

    void set_eb_velocity_bcs (amrex::Real time,
                              amrex::Vector< amrex::MultiFab* > const& vel_g_on_EB);

    void set_eb_temperature_bcs ( int const a_lev, amrex::Real const a_time,
                                  amrex::Vector<amrex::MultiFab*> a_Teb,
                                  amrex::MultiFab const* a_epf,
                                  amrex::MultiFab const* a_Tf,
                                  int const a_avgTp_comp,
                                  amrex::MultiFab const* a_avg_data);

    void set_eb_scalar_bcs ( amrex::Vector< amrex::MultiFab* > const& scalars_on_EB,
                             amrex::Vector< amrex::MultiFab* > const& species_on_EB,
                             int const a_has_enthalpy, int const a_nspecies);

    void set_gradp_bcs (const amrex::Box& bx,
                        const int lev,
                        amrex::FArrayBox& gp_fab,
                        amrex::Box& domain);


    void calc_bc_areas (int const lev,
                        amrex::BoxArray            const& a_grids,
                        amrex::DistributionMapping const& a_dmap,
                        amrex::EBFArrayBoxFactory* factory);

    amrex::Real calc_bc_area (int             const       a_lev,
                              int             const       a_dir,
                              int             const       a_bcv,
                              amrex::Box      const&      a_bx,
                              BC_t            const&      a_bc,
                              amrex::MultiFab const&      a_MF,
                              amrex::EBFArrayBoxFactory*  a_factory);

    amrex::Real calc_eb_bc_area (int             const       a_lev,
                                 int             const       a_bcv,
                                 BC_t            const&      a_bc,
                                 amrex::MultiFab const&      a_MF,
                                 amrex::EBFArrayBoxFactory*  a_factory);

    template <typename F1>
    AMREX_FORCE_INLINE
    void set_bcs (int const a_lev, amrex::Real /*time*/,
                  F1 bcs_function, amrex::MultiFab* mf_in)
    {
      const int ncomp = mf_in->nComp();

      using BctArr = amrex::Array4<const int>;

      BctArr const& bct_ilo = m_bc_list.bc_ilo[a_lev]->const_array();
      BctArr const& bct_jlo = m_bc_list.bc_jlo[a_lev]->const_array();
      BctArr const& bct_klo = m_bc_list.bc_klo[a_lev]->const_array();

      amrex::Array<BctArr const*, 3> bct_lo_arr = {&bct_ilo, &bct_jlo, &bct_klo};

      BctArr const& bct_ihi = m_bc_list.bc_ihi[a_lev]->const_array();
      BctArr const& bct_jhi = m_bc_list.bc_jhi[a_lev]->const_array();
      BctArr const& bct_khi = m_bc_list.bc_khi[a_lev]->const_array();

      amrex::Array<BctArr const*, 3> bct_hi_arr = {&bct_ihi, &bct_jhi, &bct_khi};

      amrex::Array<amrex::Array<BctArr const*, 3>*, 2> bct_vect = {&bct_lo_arr, &bct_hi_arr};

      amrex::Box domain(m_geom[a_lev].Domain());

      amrex::IntVect dom_lo(domain.loVect());
      amrex::IntVect dom_hi(domain.hiVect());

      amrex::Array<amrex::IntVect*, 2> dom_vect = {&dom_lo, &dom_hi};

      for (amrex::MFIter mfi(*mf_in, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        amrex::FArrayBox& mf_fab = (*mf_in)[mfi];

        const amrex::EBFArrayBox& mf_EB_fab = static_cast<amrex::EBFArrayBox const&>(mf_fab);
        const amrex::EBCellFlagFab& flags_fab = mf_EB_fab.getEBCellFlagFab();

        amrex::IntVect mf_lo(mf_fab.loVect());
        amrex::IntVect mf_hi(mf_fab.hiVect());

        amrex::Array<amrex::IntVect*, 2> mf_vect = {&mf_lo, &mf_hi};

        amrex::Array4<amrex::Real> const& mf_arr = mf_fab.array();
        amrex::Array4<const amrex::EBCellFlag> const& flags_arr = flags_fab.const_array();

        for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

          amrex::IntVect bx_lo_hi(mf_hi);
          amrex::IntVect bx_hi_lo(mf_lo);

          bx_lo_hi[dir] = dom_lo[dir]-1;
          bx_hi_lo[dir] = dom_hi[dir]+1;

          amrex::Box bx_lo(mf_lo, bx_lo_hi);
          amrex::Box bx_hi(bx_hi_lo, mf_hi);

          amrex::Array<amrex::Box*, 2> bx_vect = {&bx_lo, &bx_hi};

          for (int face(-1); face <= 1; face += 2) {

            const int face_idx = (face+1) / 2;

            const int face_flag = amrex::max(0,
                face*((*(mf_vect[face_idx]))[dir]-((*(dom_vect[face_idx]))[dir])));

            if (face_flag > 0) {

              amrex::Box bx = *(bx_vect[face_idx]);

              amrex::Array4<const int> const& bct_arr = *((*(bct_vect[face_idx]))[dir]);
              amrex::IntVect dom = *(dom_vect[face_idx]);
              amrex::ParallelFor(bx, ncomp, [bct_arr,dom,mf_arr,flags_arr,
                  bcs_function,dir,face]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                amrex::GpuArray<int, AMREX_SPACEDIM> ijk = {i,j,k};
                amrex::IntVect dom_ijk(ijk.data());
                dom_ijk[dir] = dom[dir]+face;

                const int bct = bct_arr(dom_ijk,0);
                const int bcv = bct_arr(dom_ijk,1);

                dom_ijk[dir] = dom[dir];

                bcs_function(bct, bcv, {i,j,k}, dom_ijk, n, mf_arr, flags_arr, dir);
              });
            }
          }
        }

      } // end MFIter loop

      mf_in->FillBoundary(m_geom[a_lev].periodicity());

      amrex::EB_set_covered(*mf_in, 0, mf_in->nComp(),
        mf_in->nGrow(), m_covered_val);

    } // end set_bc

    template <typename F1>
    AMREX_FORCE_INLINE
    void set_bcs_2D (int const a_lev, amrex::Real /*time*/,
                     F1 bcs_function, amrex::MultiFab* mf_in)
    {
      const int ncomp = mf_in->nComp();

      using BctArr = amrex::Array4<const int>;

      BctArr const& bct_ilo = m_bc_list.bc_ilo[a_lev]->const_array();
      BctArr const& bct_jlo = m_bc_list.bc_jlo[a_lev]->const_array();
      BctArr const& bct_klo = m_bc_list.bc_klo[a_lev]->const_array();

      amrex::Array<BctArr const*, 3> bct_lo_arr = {&bct_ilo, &bct_jlo, &bct_klo};

      BctArr const& bct_ihi = m_bc_list.bc_ihi[a_lev]->const_array();
      BctArr const& bct_jhi = m_bc_list.bc_jhi[a_lev]->const_array();
      BctArr const& bct_khi = m_bc_list.bc_khi[a_lev]->const_array();

      amrex::Array<BctArr const*, 3> bct_hi_arr = {&bct_ihi, &bct_jhi, &bct_khi};

      amrex::Array<amrex::Array<BctArr const*, 3>*, 2> bct_vect = {&bct_lo_arr, &bct_hi_arr};

      amrex::Box domain(m_geom[a_lev].Domain());

      amrex::IntVect dom_lo(domain.loVect());
      amrex::IntVect dom_hi(domain.hiVect());

      amrex::Array<amrex::IntVect*, 2> dom_vect = {&dom_lo, &dom_hi};

      for (amrex::MFIter mfi(*mf_in, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        amrex::FArrayBox& mf_fab = (*mf_in)[mfi];

        amrex::IntVect mf_lo(mf_fab.loVect());
        amrex::IntVect mf_hi(mf_fab.hiVect());

        amrex::Array<amrex::IntVect*, 2> mf_vect = {&mf_lo, &mf_hi};

        amrex::Array4<amrex::Real> const& mf_arr = mf_fab.array();

        for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

          amrex::IntVect bx_lo_lo(mf_lo), bx_lo_hi(mf_hi);
          amrex::IntVect bx_hi_lo(mf_lo), bx_hi_hi(mf_hi);

          bx_lo_lo[dir] = dom_lo[dir]-1;
          bx_lo_hi[dir] = dom_lo[dir]-1;
          bx_hi_lo[dir] = dom_hi[dir]+1;
          bx_hi_hi[dir] = dom_hi[dir]+1;

          amrex::Box bx_lo(bx_lo_lo, bx_lo_hi);
          amrex::Box bx_hi(bx_hi_lo, bx_hi_hi);

          amrex::Array<amrex::Box*, 2> bx_vect = {&bx_lo, &bx_hi};

          for (int face(-1); face <= 1; face += 2) {

            const int face_idx = (face+1) / 2;

            const int face_flag = amrex::max(0,
                face*((*(mf_vect[face_idx]))[dir] - (*(dom_vect[face_idx]))[dir]));

            if (face_flag > 0)
            {
              amrex::Box bx = *(bx_vect[face_idx]);

              amrex::Array4<const int> const& bct_arr = *((*(bct_vect[face_idx]))[dir]);
              amrex::IntVect dom = *(dom_vect[face_idx]);

              amrex::ParallelFor(bx, ncomp, [bct_arr,dom,mf_arr,dir,face,
                  bcs_function]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                amrex::GpuArray<int, AMREX_SPACEDIM> ijk = {i,j,k};
                amrex::IntVect dom_ijk(ijk.data());
                dom_ijk[dir] = dom[dir]+face;

                const int bct = bct_arr(dom_ijk,0);

                dom_ijk[dir] = dom[dir];

                amrex::IntVect near_ijk(ijk.data());
                near_ijk[dir] -= face;

                bcs_function(bct, {i,j,k}, dom_ijk, near_ijk, n, mf_arr);
              });
            }
          }
        }

      } // end MFIter loop

      mf_in->FillBoundary(m_geom[a_lev].periodicity());

      amrex::EB_set_covered(*mf_in, 0, mf_in->nComp(),
          mf_in->nGrow(), m_covered_val);
    }

    const amrex::Vector<int>& bc_xlo () const { return m_bc_lo[0]; }
    const amrex::Vector<int>& bc_xhi () const { return m_bc_hi[0]; }
    const amrex::Vector<int>& bc_ylo () const { return m_bc_lo[1]; }
    const amrex::Vector<int>& bc_yhi () const { return m_bc_hi[1]; }
    const amrex::Vector<int>& bc_zlo () const { return m_bc_lo[2]; }
    const amrex::Vector<int>& bc_zhi () const { return m_bc_hi[2]; }

    int bc_xlo (const int i) const { return m_bc_lo[0][i]; }
    int bc_xhi (const int i) const { return m_bc_hi[0][i]; }
    int bc_ylo (const int i) const { return m_bc_lo[1][i]; }
    int bc_yhi (const int i) const { return m_bc_hi[1][i]; }
    int bc_zlo (const int i) const { return m_bc_lo[2][i]; }
    int bc_zhi (const int i) const { return m_bc_hi[2][i]; }

    // Flag to solve fluid equations
    int delp_dir () const { return m_delp_dir; }
    const amrex::RealVect& delp () const { return m_delp; }
    amrex::Real delp (const int i) const { return m_delp[i]; }

    void set_delp (const int i, const amrex::Real val) { m_delp[i] = val; }

    // Flags for DEM domain extent reinforcement.
    const std::array<int,6>& domain_bc () const { return m_domain_bc; }
    int domain_bc (const int i) const { return m_domain_bc[i]; }

    const std::bitset<6>& ls_wall_flag () const { return m_ls_wall_flag; }

    // Use AMReX LinOpBCType to define the domain extents
    const std::array<amrex::LinOpBCType,3>& ppe_lobc () const { return m_ppe_lobc; }
    const std::array<amrex::LinOpBCType,3>& ppe_hibc () const { return m_ppe_hibc; }

    const std::array<amrex::LinOpBCType,3>& diff_vel_lobc () const { return m_diff_vel_lobc; }
    const std::array<amrex::LinOpBCType,3>& diff_vel_hibc () const { return m_diff_vel_hibc; }

    const std::array<amrex::LinOpBCType,3>& diff_scal_lobc () const { return m_diff_scal_lobc; }
    const std::array<amrex::LinOpBCType,3>& diff_scal_hibc () const { return m_diff_scal_hibc; }

    const std::array<amrex::LinOpBCType,3>& diff_temperature_lobc () const { return m_diff_temperature_lobc; }
    const std::array<amrex::LinOpBCType,3>& diff_temperature_hibc () const { return m_diff_temperature_hibc; }

    const std::array<amrex::LinOpBCType,3>& diff_species_lobc () const { return m_diff_species_lobc; }
    const std::array<amrex::LinOpBCType,3>& diff_species_hibc () const { return m_diff_species_hibc; }

    const std::array<amrex::LinOpBCType,3>& mac_lobc () const { return m_mac_lobc; }
    const std::array<amrex::LinOpBCType,3>& mac_hibc () const { return m_mac_hibc; }

    // Data structure for BC data
    const amrex::Vector<BC_t>& bc () const { return m_bc; }
    BC_t& bc (const int bcv) { return m_bc[bcv]; }
    const BC_t& bc (const int bcv) const { return m_bc[bcv]; }

    /// Number of boundary conditions
    int count () const noexcept { return m_bc.size(); }

    amrex::Real& h_bc_u_g (const int bcv) { return m_h_bc_u_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_u_g () { return m_h_bc_u_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_u_g () { return m_bc_u_g; }

    amrex::Real& h_bc_v_g (const int bcv) { return m_h_bc_v_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_v_g () { return m_h_bc_v_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_v_g () { return m_bc_v_g; }

    amrex::Real& h_bc_w_g (const int bcv) { return m_h_bc_w_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_w_g () { return m_h_bc_w_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_w_g () { return m_bc_w_g; }

    amrex::Real& h_bc_T (const int bcv) { return m_h_bc_T[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_T () { return m_h_bc_T; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_T () { return m_bc_T; }

    amrex::Real& h_bc_h (const int bcv) { return m_h_bc_h[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_h () { return m_h_bc_h; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_h () { return m_bc_h; }

    amrex::Real& h_bc_rho (const int bcv) { return m_h_bc_rho[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_rho () { return m_h_bc_rho; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_rho () { return m_bc_rho; }

    amrex::Real& h_bc_tracer (const int bcv) { return m_h_bc_tracer[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_tracer () { return m_h_bc_tracer; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_tracer () { return m_bc_tracer; }

    amrex::Real& h_bc_epf (const int bcv) { return m_h_bc_epf[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_epf () { return m_h_bc_epf; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_epf () { return m_bc_epf; }

    amrex::Real& h_bc_p (const int bcv) { return m_h_bc_p[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_p () { return m_h_bc_p; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_p () { return m_bc_p; }

    amrex::Real& h_bc_Xk (const int bcv) { return m_h_bc_Xk[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_Xk () { return m_h_bc_Xk; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_Xk () { return m_bc_Xk; }


    amrex::Real& get_bc_area (const int bcv);

    amrex::Real& get_bc_fab_area (const int a_bcv, std::pair<int,int> a_index);

    int& get_dir (const int bcv) { return m_dir[bcv]; }

    BCList const& get_bc_list () const { return m_bc_list; }

    amrex::Vector<amrex::BCRec> const& get_velocity_bcrec () const noexcept
    { return m_bcrec_velocity; }

    amrex::Vector<amrex::BCRec> const& get_hydro_velocity_bcrec () const noexcept
    { return m_bcrec_hydro_velocity; }

    amrex::Vector<amrex::BCRec> const& get_volfrac_bcrec () const noexcept
    { return m_bcrec_volfrac; }

    amrex::Vector<amrex::BCRec> const& get_density_bcrec () const noexcept
    { return m_bcrec_density; }

    amrex::Vector<amrex::BCRec> const& get_enthalpy_bcrec () const noexcept
    { return m_bcrec_enthalpy; }

    amrex::Vector<amrex::BCRec> const& get_tracer_bcrec () const noexcept
    { return m_bcrec_tracer; }

    amrex::Vector<amrex::BCRec> const& get_species_bcrec () const noexcept
    { return m_bcrec_species; }

    // The GPU  versions of BCRec are needed for calls to AMReX::Hydro
    amrex::BCRec const* get_velocity_bcrec_device_ptr () const noexcept
    { return m_bcrec_velocity_d.data(); }

    amrex::BCRec const* get_hydro_velocity_bcrec_device_ptr () const noexcept
    { return m_bcrec_hydro_velocity_d.data(); }

    amrex::BCRec const* get_volfrac_bcrec_device_ptr () const noexcept
    { return m_bcrec_volfrac_d.data(); }

    amrex::BCRec const* get_density_bcrec_device_ptr () const noexcept
    { return m_bcrec_density_d.data(); }

    amrex::BCRec const* get_enthalpy_bcrec_device_ptr () const noexcept
    { return m_bcrec_enthalpy_d.data(); }

    amrex::BCRec const* get_tracer_bcrec_device_ptr () const noexcept
    { return m_bcrec_tracer_d.data(); }

    amrex::BCRec const* get_species_bcrec_device_ptr () const noexcept
    { return m_bcrec_species_d.data(); }

    const amrex::Geometry& Geom (int a_lev) noexcept { return m_geom[a_lev]; }

    void fillpatch ( int const a_lev, amrex::Real const a_time, BCFillVar a_var,
                     amrex::Vector<amrex::MultiFab*> const& a_MF,
                     int const a_ngrow );

    void fillpatch ( int const a_lev, amrex::Real const a_time, BCFillVar a_var,
                     amrex::Vector<amrex::MultiFab*> const& a_fill_data,
                     amrex::MultiFab* a_lev_MF, int const a_ngrow );

    void fillcoarsepatch ( int const a_lev, amrex::Real const a_time, BCFillVar a_var,
                           amrex::Vector<amrex::MultiFab*> const& a_fill_data,
                           amrex::MultiFab* a_lev_MF, int const a_ngrow );

    void get_fill_bcs ( BCFillVar a_var, int const a_ncomp,
                        amrex::Vector<amrex::BCRec>& a_bcrec,
                        amrex::Real*& a_bc_values );

    MFIXEmbeddedBoundaries& eb_parms () { return m_embedded_boundaries; }

  private:

    int m_verbose = 0;

    int const m_max_level;

    amrex::Vector<amrex::Geometry> const& m_geom;

    amrex::Vector<amrex::IntVect> m_ref_ratio;

    BCList& m_bc_list;

    // BC types for each different variable
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_vel_bc_types;
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_rho_bc_types;
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_T_bc_types;
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_trac_bc_types;
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_Xk_bc_types;

    // BCs
    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_u_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_u_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_v_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_v_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_w_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_w_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_T;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_T;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_h;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_h;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_rho;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_rho;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_tracer;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_tracer;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_epf;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_epf;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_p;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_p;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_Xk;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_Xk;

    // Lists of BCs applied to the domain extent
    std::array<amrex::Vector<int>,3> m_bc_lo;
    std::array<amrex::Vector<int>,3> m_bc_hi;

    // Domain face of BC
    // lo x/y/z: 0, 2, 4
    // hi x/y/z: 1: 3, 5
    amrex::Vector<int> m_dir;

    // Flag to solve fluid equations
    int m_delp_dir = -1;
    amrex::RealVect m_delp;

    // Flags for DEM domain extent reinforcement.
    std::array<int,6> m_domain_bc;

    std::bitset<6> m_ls_wall_flag;

    // Use AMReX LinOpBCType to define the domain extents
    std::array<amrex::LinOpBCType,3> m_ppe_lobc;
    std::array<amrex::LinOpBCType,3> m_ppe_hibc;

    std::array<amrex::LinOpBCType,3> m_diff_vel_lobc;
    std::array<amrex::LinOpBCType,3> m_diff_vel_hibc;

    std::array<amrex::LinOpBCType,3> m_diff_scal_lobc;
    std::array<amrex::LinOpBCType,3> m_diff_scal_hibc;

    std::array<amrex::LinOpBCType,3> m_diff_temperature_lobc;
    std::array<amrex::LinOpBCType,3> m_diff_temperature_hibc;

    std::array<amrex::LinOpBCType,3> m_diff_species_lobc;
    std::array<amrex::LinOpBCType,3> m_diff_species_hibc;

    std::array<amrex::LinOpBCType,3> m_mac_lobc;
    std::array<amrex::LinOpBCType,3> m_mac_hibc;

    MFIXEmbeddedBoundaries& m_embedded_boundaries;

    amrex::Vector<amrex::Real> m_area;

    amrex::Vector<std::map<std::pair<int,int>, amrex::Real>> m_fab_area;

    // Data structure for BC data
    amrex::Vector<BC_t> m_bc;

    amrex::Real const m_covered_val = 1.e+40;

    amrex::Vector<amrex::BCRec> m_bcrec_velocity;
    amrex::Vector<amrex::BCRec> m_bcrec_volfrac;
    amrex::Vector<amrex::BCRec> m_bcrec_density;
    amrex::Vector<amrex::BCRec> m_bcrec_enthalpy;
    amrex::Vector<amrex::BCRec> m_bcrec_tracer;
    amrex::Vector<amrex::BCRec> m_bcrec_species;

    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_velocity_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_volfrac_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_density_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_enthalpy_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_tracer_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_species_d;

    // This is a special set of bc's that replace foextrap by extdir
    // for *normal* velocity only at faces that are identified as
    // pressure inflow (pinf_).   We need this because otherwise the
    // hydro routines will enforce "no backflow at outflow" and will
    // zero the inflow velocity at a pressure inflow
    amrex::Vector<amrex::BCRec> m_bcrec_hydro_velocity;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_hydro_velocity_d;
};


struct BC_t {

  // Constructor
  BC_t (MFIXFluidPhase* fluid_in) :
    fluid(fluid_in)
  {}

  // BC type (mi, pi, po, eb)
  int type;

  const amrex::RealBox* region;

  int Lo ( int const a_dir, amrex::Real a_plo, amrex::Real a_dx ) const noexcept
  { return static_cast<int>(amrex::Math::floor((region->lo(a_dir)-a_plo)/a_dx + 0.5)); }

  int Hi ( int const a_dir, amrex::Real a_plo, amrex::Real a_dx ) const noexcept
  { return static_cast<int>(amrex::Math::floor((region->hi(a_dir)-a_plo)/a_dx + 0.5))-1; }

  EB_t eb;

  FLUID_t fluid;

  amrex::Vector<SOLIDS_t> solids;

  const SOLIDS_t* get_solid (const int phase) const
  {
    for (int n(0); n < solids.size(); ++n)
      if (phase == solids[n].phase)
        return &(solids[n]);

    amrex::Abort("Error: solid type not found");
    return nullptr;
  }
};

#endif
