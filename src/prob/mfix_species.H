#ifndef _MFIX_SPECIES_H_
#define _MFIX_SPECIES_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>


namespace species_algorithms {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int linear_search (const char* names_array,
                   const int* names_idxs,
                   const int size,
                   const char* name)
{
  int name_size(0);
  while (name[name_size] != '\0')
  {
    AMREX_ASSERT(name_size < 80);
    name_size++;
  }
  AMREX_ASSERT(name_size > 0);

  for (int i(0); i < size; ++i) {

    const char* name_string = &(names_array[names_idxs[i]]);
    const int length = names_idxs[i+1] - names_idxs[i];

    if (length == name_size) {
      int match = 1;
      for (int n(0); n < length; ++n) {
        if (name_string[n] != name[n]) {
          match = 0;
          break;
        }
      }

      if (match)
        return i;
    }
  }

  return -1;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int binary_search (const char* names_array,
                   const int* names_idxs,
                   const int size,
                   const char* name)
{
  int left = 0;
  int right = size-1;

  int name_size(0);
  while (name[name_size] != '\0')
  {
    AMREX_ASSERT(name_size < 80);
    name_size++;
  }
  AMREX_ASSERT(name_size > 0);

  while (left <= right) {
    const int mid = left + (right - left) / 2;

    const char* name_string = &(names_array[names_idxs[mid]]);
    const int length = names_idxs[mid+1] - names_idxs[mid];

    if (length == name_size) {
      int match = 1;
      for (int n(0); n < length; ++n) {
        if (name_string[n] != name[n]) {
          match = 0;
          break;
        }
      }

      if (match)
        return mid;
    }

    {
      for (int n(0); n < amrex::max(length, name_size); ++n) {
        if (n >= length) {
          left = mid+1;
          break;
        } else if (n >= name_size) {
          right = mid-1;
          break;
        }

        if (name_string[n] == name[n]) {
          continue;
        } else if (name_string[n] < name[n]) {
          left = mid+1;
          break;
        } else {
          right = mid-1;
          break;
        }
      }
    }
  }

  return -1;
}

} // end namespace species_algorithms


//// Forward declaration
//class MFIXSpeciesParms;

class MFIXSpecies
{
  public:
    //
    struct DiffusivityModel { enum {Undefined=-1, Constant}; };

    // Constructor
    MFIXSpecies();

    void Initialize ();

    // List of recognized species diffusivity models
    int diffusivity_model () const { return m_diffusivity_model; }

    // Flag to solve species equations
    int solve () const { return m_solve; }

    // Total number of species
    int nspecies () const { return m_nspecies; }

    // Specified species names
    const amrex::Vector<std::string>& names () const { return m_names; }
    const std::string& names (const int k) const { return m_names[k]; }

    // Species unique identifying code (at the moment = their index in the input
    // entries)
    const amrex::Vector<int>& IDs () const { return m_IDs; }
    int IDs (const int k) const { return m_IDs[k]; }

    // Specified species constant molecular weights
    const amrex::Vector<amrex::Real>& MW_k () const { return m_MW_k; }
    amrex::Real MW_k (const int k) const { return m_MW_k[k]; }

    // Specified species constant molecular weights
    const amrex::Vector<amrex::Real>& densities_k () const { return m_densities_k; }
    amrex::Real density_k (const int k) const { return m_densities_k[k]; }

    // Specified species constant diffusion coefficients
    amrex::Real diffusivity () const { return m_D; }

    int isInitialized () const { return m_is_initialized; }

//    const MFIXSpeciesParms& parameters () const { return *m_parameters; }

  private:
    // List of recognized species diffusivity models
    int m_diffusivity_model;

    // Flag to solve species equations
    int m_solve;

    // Total number of species
    int m_nspecies;

    // Specified species names
    amrex::Vector<std::string> m_names;

    // Species unique identifying code (at the moment = their index in the input
    // entries)
    amrex::Vector<int> m_IDs;

    // Specified species constant molecular weights
    amrex::Vector<amrex::Real> m_MW_k;

    //
    amrex::Vector<amrex::Real> m_densities_k = {};

    // Specified species constant diffusion coefficients
    amrex::Real m_D;

//    // Species parameters for GPU access
//    MFIXSpeciesParms* m_parameters;

    int m_is_initialized;
};


//class MFIXSpeciesParms
//{
//  public:
//    AMREX_GPU_HOST_DEVICE
//    MFIXSpeciesParms ();
//
//    template <amrex::RunOn run_on>
//    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//    amrex::Real MW (const int n) const
//    {
//      if (run_on == amrex::RunOn::Gpu)
//        return m_d_MW_k[n];
//      else if (run_on == amrex::RunOn::Cpu)
//        return m_h_MW_k[n];
//      else
//        amrex::Abort("Unrecognized RunOn option");
//
//      return 0;
//    }
//
//    amrex::Real* m_h_MW_k;
//    amrex::Real* m_d_MW_k;
//};


struct SPECIES_t {
  // Constructor
  SPECIES_t() :
    mass_fraction(0.)
  {}

  // Species mass fraction, used for ICs and BCs
  amrex::Real mass_fraction;
};

#endif
