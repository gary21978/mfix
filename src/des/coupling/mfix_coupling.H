#ifndef _MFIX_COUPLING_H_
#define _MFIX_COUPLING_H_

#include <AMReX_MultiFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Geometry.H>

#include <mfix_fluid.H>
#include <mfix_pc.H>

namespace GridToParticleIndexes {

/**
 * @brief Structure to hold interpolation transfer indexes.
 *
 * This structure contains indexes for physical quantities such as
 * fluid velocity and pressure gradient when they are mapped into
 * the MultiFab used to compute grid forces acting on particles.
 */
struct TxfrToParicle {
  /**
   * @brief Constructor to initialize the interpolation transfer indexes.
   *
   * @param a_include_Tf Flag to include fluid temperature.
   * @param a_include_divtau Flag to include viscous stress.
   * @param a_include_vm Flag to include fluid acceleration.
   */
  TxfrToParicle ( int const a_include_Tf
                , int const a_include_divtau
                , int const a_include_vm)
    : ncomp_Tf(a_include_Tf         ? 1 : 0)
    , ncomp_divtau(a_include_divtau ? 3 : 0)
    , ncomp_DufDt(a_include_vm      ? 3 : 0)
    // Set indices for interpolation array, velocity and pressure gradient
    // are always included. All other variables depend on model settings.
    , Tf(ncomp_vel + ncomp_grad_p)
    , divtau(Tf    + ncomp_Tf)
    , DufDt(divtau + ncomp_divtau)
    , count(DufDt  + ncomp_DufDt)
  {}

  const int ncomp_vel    = 3; ///< Number of fluid velocity components
  const int ncomp_grad_p = 3; ///< Number of pressure gradient components
  const int ncomp_Tf;         ///< Number of fluid temperature components
  const int ncomp_divtau;     ///< Number of viscous stress components
  const int ncomp_DufDt;      ///< Number of fluid acceleration components

  const int vel    =  0; ///< Index for fluid velocity.
  const int grad_p =  3; ///< Index for pressure gradient.
  const int Tf;          ///< Index for fluid temperature.
  const int divtau;      ///< Index for viscous stress.
  const int DufDt;       ///< Index for fluid acceleration.

  const int count;       ///< Total count of included quantities.

  /// Identify if fluid temperature is included: 1 yes; 0 No
  int include_Tf () const noexcept { return (ncomp_Tf > 0) ? 1 : 0; }

  /// Identify if viscous stress is included: 1 yes; 0 No
  int include_divtau () const noexcept { return (ncomp_divtau > 0) ? 1 : 0; }

  /// Identify if fluid acceleration is included: 1 yes; 0 No
  int include_vm () const noexcept { return (ncomp_DufDt > 0) ? 1 : 0; }
};

} // end namespace GridToParticleIndexes


// Drag type enumeration
struct DragModel {

  public:

    void setWenYu () { m_drag_model = m_WenYu; }
    void setGidaspow () { m_drag_model = m_Gidaspow; }
    void setBVK2 () { m_drag_model = m_BVK2; }
    void setSyamOBrien() { m_drag_model = m_SyamOBrien; }
    void setUserDrag() { m_drag_model = m_UserDrag; }

    void setSyamOBrien_coeffs ( amrex::Real const a_c1,
                                amrex::Real const a_d1 ) {
      m_SyamOBrien_coeff_c1 = a_c1;
      m_SyamOBrien_coeff_d1 = a_d1;
    }

    int WenYu () const { return (m_drag_model == m_WenYu ? 1 : 0); }
    int Gidaspow () const { return (m_drag_model == m_Gidaspow ? 1 : 0); }
    int BVK2 () const { return (m_drag_model == m_BVK2 ? 1 : 0); }
    int SyamOBrien() const { return (m_drag_model == m_SyamOBrien ? 1 : 0); }
    int UserDrag() const { return m_drag_model == m_UserDrag ? 1 : 0; }

    amrex::Real SyamOBrien_coeff_c1() const { return m_SyamOBrien_coeff_c1; }
    amrex::Real SyamOBrien_coeff_d1() const { return m_SyamOBrien_coeff_d1; }

  private:

    enum {m_Undefined, m_WenYu, m_Gidaspow, m_BVK2, m_SyamOBrien, m_UserDrag};

    int m_drag_model = m_Undefined;

    amrex::Real m_SyamOBrien_coeff_c1 = -1.;
    amrex::Real m_SyamOBrien_coeff_d1 = -1.;
};


struct VirtualMassModel {

    void setNull () { m_virtual_mass_model = m_Null; }
    void setZuber () { m_virtual_mass_model = m_Zuber; }
    void setNijssen () { m_virtual_mass_model = m_Nijssen; }

    void setConstant ( amrex::Real a_const_coeff ) {
      m_virtual_mass_model = m_Constant;
      m_const_coeff = a_const_coeff;
    }

    int Null () const { return m_virtual_mass_model == m_Null ? 1 : 0; }
    int Constant () const { return m_virtual_mass_model == m_Constant ? 1 : 0; }
    int Zuber () const { return m_virtual_mass_model == m_Zuber ? 1 : 0; }
    int Nijssen () const { return m_virtual_mass_model == m_Nijssen ? 1 : 0; }

    /// return true if virtual mass force is included in fluid-solids coupling
    int enabled () const { return (m_virtual_mass_model == m_Undefined) ? 0 : 1; }

    amrex::Real getConstantCoeff() const { return m_const_coeff; }

  private:

    enum {m_Undefined, m_Null, m_Constant, m_Zuber, m_Nijssen};

    int m_virtual_mass_model = m_Undefined;

    amrex::Real m_const_coeff = 0.0;
};

// Convection type enumeration
struct ConvectionModel {

  public:

    void setRanzMarshall () { m_convection_model = m_RanzMarshall; }
    void setGunn () { m_convection_model = m_Gunn; }
    void setNullConvection () { m_convection_model = m_NullConvection; }

    int RanzMarshall () const { return (m_convection_model == m_RanzMarshall ? 1 : 0); }
    int Gunn () const { return (m_convection_model == m_Gunn ? 1 : 0); }
    int NullConvection () const { return (m_convection_model == m_NullConvection ? 1 : 0); }

  private:

    enum {m_Undefined, m_RanzMarshall, m_Gunn, m_NullConvection};

    int m_convection_model = m_Undefined;

};


// Convection type enumeration
struct DepositionFilter {

  public:

    void setConstantSize ( amrex::Real a_filter_size) {
      m_filter_size = a_filter_size;
      m_filter_type = m_ConstantSize;
     }

    void setVariableSize ( amrex::Real a_filter_size) {
      m_filter_size = a_filter_size;
      m_filter_type = m_VariableSize;
     }

    void setMinSolidsVolfrac( amrex::Real const a_mineps)
    { m_filter_mineps = a_mineps; }

    int enabled () const { return m_filter_type != m_None; }

    amrex::Real size () const { return m_filter_size; }
    amrex::Real mineps () const { return m_filter_mineps; }

    int ConstantSize () const { return m_filter_type == m_ConstantSize; }
    int VariableSize () const { return m_filter_type == m_VariableSize; }

  private:

    enum {m_None, m_ConstantSize, m_VariableSize };

    int m_filter_type = m_None;

    amrex::Real m_filter_size = std::numeric_limits<amrex::Real>::max();
    amrex::Real m_filter_mineps = 1.0e-6;
};


/**
 * \class CouplingOps
 * \brief Options that control fluid-solids coupling
 */
class CouplingOp
{

  public:

    CouplingOp () = default;
    ~CouplingOp () { if (m_interp_idxs != nullptr) delete m_interp_idxs; }

    /**
    * Process ParmParse variables, perform basic checks on provided settings,
    * and set runtime flags.
    */
    int Initialize ( int const a_couple_energy,
                     int const a_nspecies,
                     int const a_is_mixture,
                     int const a_solve_fluid_and_solids );

    /// Provides access to drag model settings.
    DragModel const* getDragModel () { return &m_drag_model; }

    /// Provides access to virtual mass model settings.
    VirtualMassModel const* getVirtualMassModel () { return &m_virtual_mass_model; }

    /// Provides access to convection model settings.
    ConvectionModel const* getConvectionModel () { return &m_convection_model; }

    /// Provides access to deposition filter settings.
    DepositionFilter const* getDepositionFilter () { return &m_depop_filter; }

    /// returns true if deposition filtering is enabled
    int FilterDeposition () const { return m_depop_filter.enabled(); }

    void set_include_divtau ( int const a_include_divtau )
    { m_include_divtau = a_include_divtau; }

    /// return true if viscous stress is included in fluid-solids coupling
    int include_divtau () const { return m_include_divtau; }

    /// return true if virtual mass is included in fluid-solids coupling
    int include_virtual_mass () const { return m_virtual_mass_model.enabled(); }

    amrex::Real small_number () const { return m_small_number; }
    amrex::Real large_number () const { return m_large_number; }

    int interp_ng ( ) const { return 1; }

    int interp_comps ( ) const {
      AMREX_ASSERT( m_interp_idxs != nullptr );
      return m_interp_idxs->count;
    }

    int include_energy ( ) const {
      AMREX_ASSERT( m_interp_idxs != nullptr );
      return m_interp_idxs->include_energy;
    }

    int include_species ( ) const {
      AMREX_ASSERT( m_interp_idxs != nullptr );
      return ( include_energy() && m_interp_idxs->is_mixture );
    }

    int nspecies ( ) const {
      AMREX_ASSERT( m_interp_idxs != nullptr );
      return m_interp_idxs->nspecies;
    }

    int OnSameGrids ( int const a_lev ) const {
      AMREX_ASSERT( a_lev >= 0 && a_lev < m_on_same_grids.size() );
      return m_on_same_grids[a_lev];
    }

    amrex::Vector< amrex::MultiFab* > get_interp () const noexcept {
      amrex::Vector< amrex::MultiFab*> r;
      r.reserve(m_interp.size());
      for (auto& lev_interp : m_interp) { r.push_back( lev_interp ); }
      return r;
    }

    void setup ( amrex::Vector<amrex::Geometry>& a_geom,
                 amrex::Vector< amrex::MultiFab* > const& a_epf,
                 amrex::Vector< amrex::MultiFab* > const& a_rho,
                 amrex::Vector< amrex::MultiFab* > const& a_vel,
                 amrex::Vector< amrex::MultiFab* > const& a_Tf,
                 amrex::Vector< amrex::MultiFab* > const& a_hf,
                 amrex::Vector< amrex::MultiFab* > const& a_Xfk,
                 amrex::Vector<amrex::EBFArrayBoxFactory const*> const& a_particle_factory);

    void calc_transfer_coeffs ( MFIXFluidPhase const& a_fluid,
                                MFIXParticleContainer* a_pc,
                                amrex::GeometryData const& a_geomData);

    template <typename F1>
    void calc_transfer_coeffs ( MFIXFluidPhase const& a_fluid,
                                MFIXParticleContainer* a_pc,
                                amrex::GeometryData const& a_geomData,
                                F1 DragCoeff);

    template <typename F1, typename F2>
    void calc_transfer_coeffs ( MFIXFluidPhase const& a_fluid,
                                MFIXParticleContainer* a_pc,
                                amrex::GeometryData const& a_geomData,
                                F1 DragCoeff, F2 VirtualMassCoeff);

    template <typename F1, typename F2, typename F3>
    void calc_transfer_coeffs ( MFIXFluidPhase const& a_fluid,
                                MFIXParticleContainer* a_pc,
                                amrex::GeometryData const& a_geomData,
                                F1 DragCoeff, F2 VirtualMassCoeff, F3 ConvectionCoeff);

    int reset ();

  private:

    DragModel m_drag_model;
    VirtualMassModel m_virtual_mass_model;
    ConvectionModel m_convection_model;
    DepositionFilter m_depop_filter;

    /// Flag to include shear stress in fluid-particle coupling
    int m_include_divtau = 0;

    amrex::Real const m_small_number = 1.0e-15;
    amrex::Real const m_large_number = 1.0e+32;

    DualGridAuxIndexes* m_interp_idxs;

    amrex::Vector< int > m_on_same_grids;

    int const m_interp_ng = 1;

    // Auxiliary MultiFabs
    amrex::Vector< amrex::MultiFab* > m_interp;

};
#endif
