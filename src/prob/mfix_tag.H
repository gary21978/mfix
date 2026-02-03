#ifndef _MFIX_TAG_H_
#define _MFIX_TAG_H_

#include <AMReX_GpuContainers.H>
#include <mfix_regions.H>

template <typename T>
struct MFIXTagType {

  bool operator() () const { return m_tagging; }

  bool less_than () const { return m_less_than; }

  T value ( int a_lev ) const { return m_values.at(a_lev); }

  void set ( bool a_flag, bool a_less_than,
             amrex::Vector<T> a_values ) {
    m_tagging = a_flag;
    m_less_than = a_less_than;
    m_values = a_values;
  }

  private:

  /// Flag that we are tagging this variable
  bool m_tagging = false;

  /// Use less than operator and not greater than
  bool m_less_than = false;

  /// Values used for tagging
  amrex::Vector<T> m_values;
};

class MFIXTag {

  public:

    MFIXTag () = default;
    ~MFIXTag () = default;

    void Initialize ( int const a_max_level,
                      int const a_solve_fluid,
                      int const a_solve_solids,
                      const MFIXRegions& a_regions );

    MFIXTagType<amrex::Real> rho;
    MFIXTagType<amrex::Real> grad_rho;

    MFIXTagType<amrex::Real> vel;
    MFIXTagType<amrex::Real> grad_vel;

    MFIXTagType<amrex::Real> vorticity;

    amrex::Gpu::HostVector<amrex::RealBox> h_regions;

    int num_regions () const { return m_num_regions; }

  private:

    int m_tags = 0;

    int m_num_regions = 0;
};


#endif

