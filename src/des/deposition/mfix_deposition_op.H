#ifndef MFIX_DEPOSITION_OP_H_
#define MFIX_DEPOSITION_OP_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBCellFlag.H>

#include <mfix_bc.H>
#include <mfix_pc.H>

class MFIXDepositionOp
{
  public:

    using ParticleType = MFIXParticleContainer::ParticleType;
    using PairIndex = MFIXParticleContainer::PairIndex;
    using MFIXParIter = MFIXParticleContainer::MFIXParIter;

    MFIXDepositionOp ( BCList const& a_bc_list,
                       MFIXParticleContainer* a_pc)
    : m_bc_list(a_bc_list)
    , m_pc(a_pc)
    { }

    MFIXDepositionOp ( BCList const& a_bc_list,
                       MFIXParticleContainer* a_pc,
                       const amrex::EBFArrayBoxFactory* a_lev0_factory,
                       amrex::Vector <amrex::MultiFab*> a_dst );

    // Destructor
    virtual ~MFIXDepositionOp ();

    virtual void Deposit ( )
    { amrex::Abort("MFIXDepositionOp base class has no deposit() operation!"); }

    void ApplyBC ( int const a_lev, amrex::MultiFab* a_MF);

    void Redistribute ( int const a_lev = 0 );

    void Redistribute_MaxPack ( int const a_lev );

    void Redistribute_SRD ( int const a_lev );

    void CopyToDest ( );

  protected:

    int getInputParams ();

    int OnSameGrids ( int const a_lev ) const {
      AMREX_ASSERT( a_lev >= 0 && a_lev < m_on_same_grids.size() );
      return m_on_same_grids[a_lev];
    }

    // Set option to reduce all data to level-0. This option is true
    // by default.
    void reduceToLev0 ( int const a_reduce )
    { m_reduce_to_lev0 = (a_reduce == 0 ? 0 : 1); }

    // Deposition scheme enumeration
    enum DepositionScheme {trilinear, square_dpvm, true_dpvm, centroid};

    DepositionScheme m_deposition_scheme;

    amrex::Real m_deposition_scale_factor;

    /// Solids redistribution algorithm: MaxPack or StateRedist
    enum RedistType {m_MaxPack, m_StateRedist};
    int m_redist_type;

    /**
    * When solids redistribution type is set to MaxPack, then excess solids volume
    * is redistributed to adjacent cells when the solids volume fraction of a cut-cell
    * exceeds m_redist_max_pack.
    */
    amrex::Real m_redist_max_pack;

    /**
    * When solids redistribution type is set to StateRedist, then the solids volume
    * is redistributed for cut-cells with geometric volume fractions less than
    * m_redist_vfrac.
    */
    amrex::Real m_redist_vfrac;

    int m_reduce_to_lev0 = 1;

    BCList const& m_bc_list;

    MFIXParticleContainer* m_pc;

    const amrex::EBFArrayBoxFactory* m_lev0_factory;

    amrex::Vector< int > m_on_same_grids;

    amrex::Vector< amrex::MultiFab* > m_dst;

    amrex::Vector< amrex::MultiFab* > m_ptr;
    amrex::Vector< amrex::MultiFab* > m_eps;

};


class MFIXDepOpVolAvg : public MFIXDepositionOp
{
  public:

    MFIXDepOpVolAvg ( BCList const& a_bc_list,
                      MFIXParticleContainer* a_pc,
                      const amrex::EBFArrayBoxFactory* a_lev0_factory,
                      amrex::Vector <amrex::MultiFab* > a_dst );

    ~MFIXDepOpVolAvg () = default;

    void setDepComp ( int const a_comp );

    void setDepComp ( amrex::Gpu::HostVector<int> const& a_h_depComp );

    virtual void Deposit ( );

    template <typename F>
    void Deposit ( int const a_lev, F WeightFunc);

    void Average ();

    void Average ( int const a_lev, int const a_comp, int const a_ncomp );

  private:

    int m_ncomp = 0;
    int m_deposit_only_volume = 0;

    amrex::Gpu::HostVector<int>   m_h_depComp;
    amrex::Gpu::DeviceVector<int> m_d_depComp;
};


class MFIXDepOpTxfr : public MFIXDepositionOp
{
  public:

    MFIXDepOpTxfr ( BCList const& a_bc_list,
                    MFIXParticleContainer* a_pc,
                    const amrex::EBFArrayBoxFactory* a_lev0_factory,
                    amrex::Vector <amrex::MultiFab* > a_dst,
                    amrex::Real const a_dt );

    ~MFIXDepOpTxfr () = default;

    virtual void Deposit ( );

    template <typename F>
    void Deposit (int const a_lev, F WeightFunc);

    void setVirtualMass () { m_virtual_mass = 1; }

    int includeVirtualMass () { return m_virtual_mass; }

  private:

    int m_virtual_mass = 0;

    amrex::Real const m_dt;
};


#endif
