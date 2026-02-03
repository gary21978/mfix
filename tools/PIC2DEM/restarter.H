#ifndef MFIX_RESTARTER_H_
#define MFIX_RESTARTER_H_

#include <mfix.H>
#include <mfix_pc_hex_close_pack_K.H>
#include <deposition/mfix_deposition_op.H>

struct Transfer {

  Transfer (const MFIXSolidsPhase& solids)
    : idx_eps(0)
    , idx_density(1)
    , idx_vel(2)
    , idx_temp(5)
    , idx_species(idx_temp + solids.solve_enthalpy())
    , count(idx_species + solids.nspecies()*solids.solve_species())
  {}

  const int idx_eps;
  const int idx_density;
  const int idx_vel;
  const int idx_temp;
  const int idx_species;
  const int count;
};


class MFIXDepOpPIC2DEM : public MFIXDepositionOp
{
  public:

     MFIXDepOpPIC2DEM ( BCList const& a_bc_list,
                        MFIXParticleContainer* a_pc,
                        const amrex::EBFArrayBoxFactory* a_lev0_factory,
                        amrex::Vector <amrex::MultiFab*> a_dst )
     : MFIXDepositionOp(a_bc_list, a_pc, a_lev0_factory, a_dst) { }

    virtual ~MFIXDepOpPIC2DEM ( )
    { }

    virtual void Deposit ( );

    template <typename F>
    void Deposit ( int const a_lev, F WeightFunc);

};


class MFIXRestarter
{
  public:
    MFIXRestarter (const int nlev);

    ~MFIXRestarter ();

    void allocate_avgdPIC_coarse (const mfix* mfix_coarse);

    void allocate_avgdPIC_fine (const mfix* mfix_fine);

    static void change_inputs_table ();

    void set_fine_grids_from_coarse (mfix* mfix_fine,
                                     const mfix* mfix_coarse) const;

    void deposit_PIC (const mfix* mfix_coarse);

    void set_eps_ptr (amrex::MultiFab* eps_ptr,
                      const mfix* mfix_coarse) const;

    void convert_coarse_data (const mfix* mfix_coarse,
                              mfix* mfix_fine) const;

    void get_eps_coarse (const mfix* mfix_coarse);

    void free_eps_coarse ()
    {
      if (avgd_eps_coarse != nullptr)
        delete avgd_eps_coarse;
      avgd_eps_coarse = nullptr;
    }

    void get_particles_radius (const mfix* mfix_ptr);

    void generate_particles (const amrex::Geometry& geom_coarse,
                             mfix* mfix_fine) const;

    void generate_particles (const long particles_count,
                             const amrex::Box& bx,
                             MFIXParticleContainer::ParticleTileType& particles,
                             const Hex_ClosePack* hcp_vector_ptr,
                             const Long* gen_number_ptr,
                             const Long* gen_indexes_ptr,
                             const int id,
                             const int cpu) const;

    void init_particles_data (mfix* mfix_fine) const;

    int m_refinement_ratio;

    amrex::Real m_eps_tolerance;
    amrex::Real m_eps_overflow;

    amrex::Real m_inputs_pdiameter;
    amrex::Real m_inputs_pdensity;

    int m_add_thermal_noise;

    int nlev;

    amrex::Vector<amrex::MultiFab*> avgdPIC_coarse;
    amrex::Vector<amrex::MultiFab*> avgdPIC_fine;

    amrex::MultiFab* avgd_eps_coarse;
};


class TxfrAuxiliary
{
  public:
    AMREX_GPU_HOST_DEVICE
    TxfrAuxiliary () = default;

    AMREX_GPU_HOST_DEVICE
    amrex::RealVect get_coordinates_cc (const amrex::GeometryData& geom,
                                        const amrex::IntVect& indexes) const
    {
      const amrex::Real* plo = geom.ProbLo();
      const amrex::Real* dx  = geom.CellSize();

      amrex::RealVect result(0.);

      for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
        result[dir] = (0.5 + amrex::Real(indexes[dir])) * dx[dir] + plo[dir];
      }

      return result;
    }

    AMREX_GPU_HOST_DEVICE
    amrex::IntVect get_indexes_cc (const amrex::GeometryData& geom,
                                   const amrex::RealVect& coords) const
    {
      const amrex::Real* plo = geom.ProbLo();
      const amrex::Real* dx  = geom.CellSize();

      amrex::IntVect result(0);

      for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
        result[dir] = int(amrex::Math::floor((coords[dir]-plo[dir]) / dx[dir]));
      }

      return result;
    }

    AMREX_GPU_HOST_DEVICE
    void get_weights_cc (amrex::GpuArray<amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,3>,3>& weights,
                         const amrex::Array4<const amrex::Real>& volfrac_coarse_arr,
                         const amrex::Array4<const amrex::Real>& epg_coarse_arr,
                         const amrex::Box& domain,
                         const amrex::IntVect& ijk_coarse,
                         const amrex::Array4<const amrex::EBCellFlag>& flags_coarse_arr,
                         const amrex::IntVect& ijk_fine,
                         const amrex::Array4<const amrex::EBCellFlag>& flags_fine_arr) const
    {
      if ((!flags_fine_arr(ijk_fine).isCovered()) && (flags_coarse_arr(ijk_coarse).isCovered())) {

        amrex::Real sum(0.);

        for (int ii(-1); ii <= 1; ++ii)
        for (int jj(-1); jj <= 1; ++jj)
        for (int kk(-1); kk <= 1; ++kk) {
          amrex::IntVect shift(ii,jj,kk);
          if (!flags_coarse_arr(ijk_coarse+shift).isCovered() &&
              domain.contains(ijk_coarse+shift)) {

            const amrex::Real vfrac = volfrac_coarse_arr(ijk_coarse+shift);
            const amrex::Real epg   = epg_coarse_arr(ijk_coarse+shift);

            sum += vfrac*epg;
          }
        }

        for (int ii(-1); ii <= 1; ++ii)
        for (int jj(-1); jj <= 1; ++jj)
        for (int kk(-1); kk <= 1; ++kk) {
          amrex::IntVect shift(ii,jj,kk);
          if (!flags_coarse_arr(ijk_coarse+shift).isCovered() &&
              domain.contains(ijk_coarse+shift)) {

            const amrex::Real vfrac = volfrac_coarse_arr(ijk_coarse+shift);
            const amrex::Real epg   = epg_coarse_arr(ijk_coarse+shift);

            weights[ii+1][jj+1][kk+1] = (vfrac*epg) / sum;

          } else {

            weights[ii+1][jj+1][kk+1] = 0.;
          }
        }

      } else {

        for (int ii(-1); ii <= 1; ++ii)
        for (int jj(-1); jj <= 1; ++jj)
        for (int kk(-1); kk <= 1; ++kk) {

          weights[ii+1][jj+1][kk+1] = 0.;
        }

        weights[1][1][1] = 1.;
      }
    }

    AMREX_GPU_HOST_DEVICE
    amrex::Real sum_up (const amrex::GpuArray<amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,3>,3>& weights,
                        const amrex::Box& domain,
                        const amrex::Array4<const amrex::EBCellFlag>& flags_arr,
                        const amrex::Array4<const amrex::Real>& mf_arr,
                        const amrex::IntVect& ijk,
                        const int n=0) const
    {
      amrex::Real result(0.);

      for (int ii(-1); ii <= 1; ++ii)
      for (int jj(-1); jj <= 1; ++jj)
      for (int kk(-1); kk <= 1; ++kk) {
        amrex::IntVect shift(ii,jj,kk);
        if (!flags_arr(ijk+shift).isCovered() &&
            domain.contains(ijk+shift)) {
          result += weights[ii+1][jj+1][kk+1] * mf_arr(ijk+shift, n);
        }
      }

      return result;
    }


    AMREX_GPU_HOST_DEVICE
    amrex::RealVect get_coordinates_nd (const amrex::GeometryData& geom,
                                        const amrex::IntVect& indexes) const
    {
      const amrex::Real* plo = geom.ProbLo();
      const amrex::Real* dx  = geom.CellSize();

      amrex::RealVect result(0.);

      for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
        result[dir] = amrex::Real(indexes[dir]) * dx[dir] + plo[dir];
      }

      return result;
    }

    AMREX_GPU_HOST_DEVICE
    amrex::IntVect get_indexes_nd (const amrex::GeometryData& geom,
                                   const amrex::RealVect& coords) const
    {
      const amrex::Real* plo = geom.ProbLo();
      const amrex::Real* dx  = geom.CellSize();

      amrex::IntVect result(0);

      for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
        result[dir] = int(amrex::Math::round((coords[dir]-plo[dir]) / dx[dir]));
      }

      return result;
    }

};

#endif
