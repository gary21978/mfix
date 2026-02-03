#include <mfix_mlmg_options.H>

namespace MfixUtil {

MLMGOptions::MLMGOptions(const std::string &prefix) { parse_options(prefix); }

void MLMGOptions::parse_options(const std::string& prefix)
{
   amrex::ParmParse pp(prefix);

   pp.query("mg_max_coarsening_level", max_coarsening_level);
   pp.query("mg_rtol", mg_rtol);
   pp.query("mg_atol", mg_atol);

   pp.query("bottom_solver", bottom_solver_type);
   pp.query("hypre_namespace", hypre_namespace);
   pp.query("hypre_interface", hypre_interface);
}

void MLMGOptions::apply(amrex::MLMG &mlmg) const
{
   amrex::ignore_unused(mlmg);

   if (bottom_solver_type == "hypre") {
#ifdef AMREX_USE_HYPRE
      mlmg.setHypreOptionsNamespace(hypre_namespace);
      if (hypre_interface == "ij")
         mlmg.setHypreInterface(amrex::Hypre::Interface::ij);
      else if (hypre_interface == "semi_structured")
         mlmg.setHypreInterface(amrex::Hypre::Interface::semi_structed);
      else if (hypre_interface == "structured")
         mlmg.setHypreInterface(amrex::Hypre::Interface::structed);
      else
         amrex::Abort(
               "Invalid hypre interface. Valid options: ij semi_structured "
               "structured");
#else
      amrex::Abort("mfix was not built with hypre support");
#endif
   }
}

void MLMGOptions::apply(Hydro::MacProjector &mac_proj) const
{
   apply(mac_proj.getMLMG());
}

void MLMGOptions::apply(Hydro::NodalProjector& nodal_proj) const
{
   apply(nodal_proj.getMLMG());
}

}
