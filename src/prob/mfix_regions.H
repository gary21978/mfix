#ifndef _MFIXREGIONS_H_
#define _MFIXREGIONS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

// Name and extents of regions
class MFIXRegions
{
  public:
    // Constructor
    MFIXRegions ()
      : m_is_initialized(0)
    {}

    // Initialization function -- read input parameters
    void Initialize ();

    // Returns the index of region with matching name.
    // -1 is returned if the region is not matched.
    int getIndex(const std::string& input) const;

    // Return region extents by region name
    const amrex::RealBox* getRegion(const std::string& name) const;

    int isInitialized() const
    { return m_is_initialized; }

  private:
    amrex::Vector<std::string> m_names;

    amrex::Vector<amrex::RealBox> m_extents;

    int m_is_initialized;
};

#endif
