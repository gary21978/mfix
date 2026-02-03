#ifndef MLMG_OPTIONS_H_
#define MLMG_OPTIONS_H_

#include <AMReX_ParmParse.H>
#include <AMReX_MLLinOp.H>
#include <hydro_MacProjector.H>
#include <hydro_NodalProjector.H>

namespace MfixUtil {

/* Class to parse and apply necessary MLMG options to the Nodal and Mac Projectors.
 * Currently limited to options not already parsed inside AMRex */
class MLMGOptions {
   public:
      // Constructor
      explicit MLMGOptions(const std::string &prefix);

      void apply(amrex::MLMG &mlmg) const;
      void apply(Hydro::NodalProjector &nodal_proj) const;
      void apply(Hydro::MacProjector &mac_proj) const;

      //! Relative tolerance for convergence of MLMG solvers
      amrex::Real mg_rtol{1.0e-11};

      //! Absolute tolerance for convergence checks
      amrex::Real mg_atol{1.0e-14};

      //! Max coarsening level
      int max_coarsening_level{100};

      //! Bottom solver type bicgstab, cg, hypre, etc.
      std::string bottom_solver_type{"bicgstab"};

      //! Options namespace for hypre solvers
      std::string hypre_namespace{"hypre"};

   private:
      void parse_options(const std::string &prefix);

      //! Options of hypre interface
      std::string hypre_interface{"ij"};
};

}

#endif
