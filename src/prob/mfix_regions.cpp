#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#include <AMReX_ParmParse.H>
#include <mfix_regions.H>

void
MFIXRegions::Initialize ()
{
  m_is_initialized = 1;

  // Read in the region data
  amrex::ParmParse pp("mfix");

  std::vector<std::string> regions_in;
  pp.queryarr("regions", regions_in);

  for(size_t lc=0; lc < regions_in.size(); lc++){

    amrex::RealBox new_box;

    amrex::Vector<amrex::Real> reg_lo(3);
    amrex::Vector<amrex::Real> reg_hi(3);

    {
      std::string region_prefix = "regions."+regions_in[lc];
      amrex::ParmParse ptmp(region_prefix);

      ptmp.getarr("lo", reg_lo, 0, 3);
      ptmp.getarr("hi", reg_hi, 0, 3);

      new_box.setLo(reg_lo);
      new_box.setHi(reg_hi);

      m_names.push_back(regions_in[lc]);
      m_extents.push_back(new_box);
    }
  }

} // end Initialize


/*-----------------------------------------------------------------\
|  Given 'name', return the index of the index of the region with  |
|  the matching name.                                              |
\-----------------------------------------------------------------*/
int
MFIXRegions::getIndex (const std::string& name) const
{
  for (int lc=0; lc < m_names.size(); lc++) {
    if (name.compare(m_names[lc]) == 0) {
      return lc;
    }
  }

  return -1; // Failed to find a match
}


/*-----------------------------------------------------------------\
|  Given 'name', return a pointer to the RealBox of the region     |
|  with the matching name.                                         |
\-----------------------------------------------------------------*/
const amrex::RealBox*
MFIXRegions::getRegion (const std::string& name) const
{
  for (int lc=0; lc < m_names.size(); lc++) {
    if (name.compare(m_names[lc]) == 0) {
      return &m_extents[lc];
    }
  }

  return nullptr;
}
