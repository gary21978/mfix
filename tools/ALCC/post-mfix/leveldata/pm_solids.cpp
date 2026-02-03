#include <AMReX_ParmParse.H>

#include <pm_solids.H>
#include <pm_diffusion.H>

using namespace amrex;

pm_solids::
~pm_solids ()
{
  for (auto& lev_data : m_data) { lev_data.reset( nullptr ); }
}

pm_solids::
pm_solids ( int const a_max_level,
            amrex::Vector<amrex::Geometry>            const & a_geom,
            amrex::Vector<amrex::DistributionMapping> const & a_dmap,
            amrex::Vector<amrex::BoxArray>            const & a_grids,
            pm_particles* a_particles)
  : m_max_level(a_max_level)
  , m_geom(a_geom)
  , m_dmap(a_dmap)
  , m_grids(a_grids)
{
  m_variables = {"alpha_p", "Up", "Vp", "Wp", "Up^2", "Vp^2", "Wp^2",
                 "Up*Vp","Up*Wp", "Vp*Wp"};

  // Number of 'base' variables
  int ncomps = static_cast<int>(m_variables.size());

  m_comps = ncomps + 7; // 6 Pressure + granular energy

  for (int lc(0); lc < ncomps; ++lc) {
    m_variable_map[m_variables[lc]] = lc;
  }

  m_LP_to_EL_map["alpha_p"] = "volume";
  m_LP_to_EL_map["Up"] = "velx";
  m_LP_to_EL_map["Vp"] = "vely";
  m_LP_to_EL_map["Wp"] = "velz";

  m_data.resize(get_nlev());
  for (int lev(0); lev < get_nlev(); ++lev) {

    m_data[lev].reset(new MultiFab(m_grids[lev], m_dmap[lev], m_comps, m_ngrow));
    m_data[lev]->setVal(0.0, 0, m_comps, m_ngrow);

    for (int comp(0); comp<ncomps; ++comp) {

      if ( m_variables[comp].find('*') != std::string::npos) {

          auto const pos = m_variables[comp].find('*');

          std::string const EL_var1 = m_variables[comp].substr(0, pos);
          std::string const LP_var1 = m_LP_to_EL_map[EL_var1];
          AMREX_ALWAYS_ASSERT( a_particles->contains(LP_var1) );

          std::string const EL_var2 = m_variables[comp].substr(pos+1);
          std::string const LP_var2 = m_LP_to_EL_map[EL_var2];
          AMREX_ALWAYS_ASSERT( a_particles->contains(LP_var2) );

          Print() << "  depositing " << LP_var1+"*"+LP_var2
            << " to " << EL_var1+"*"+EL_var2 << "\n";

          a_particles->deposit_mult(lev, m_data[lev].get(), comp, LP_var1, LP_var2);

      } else if ( m_variables[comp].find('^') != std::string::npos) {

          auto const pos = m_variables[comp].find('^');

          std::string const EL_var = m_variables[comp].substr(0, pos);
          std::string const LP_var = m_LP_to_EL_map[EL_var];
          AMREX_ALWAYS_ASSERT( a_particles->contains(LP_var) );

          int power = std::stoi( m_variables[comp].substr(pos+1) );

          Print() << "  depositing std::pow(" << LP_var+","+std::to_string(power)+")"
            << " to " << EL_var+"^"+std::to_string(power) << "\n";

          a_particles->deposit_pow(lev, m_data[lev].get(), comp, LP_var, 1.0, power);


      } else {

          std::string const EL_var = m_variables[comp];
          std::string const LP_var = m_LP_to_EL_map[EL_var];

          Print() << "  depositing " << LP_var << " to " << EL_var << "\n";

          if (toLower(LP_var).compare("volume") )
          { AMREX_ALWAYS_ASSERT( a_particles->contains(LP_var) ); }

          a_particles->deposit(lev, m_data[lev].get(), comp, LP_var);
      }

    }
    m_data[lev]->FillBoundary(m_geom[lev].periodicity());
  }

  int const alpha_p_comp = m_variable_map["alpha"];

  // Initialize the diffuison operator.
  pm_diffusion diffusion( m_max_level, m_geom, m_dmap, m_grids);

  diffusion.reset(GetVecOfConstPtrs(m_data), alpha_p_comp,
                  a_particles->get_mean_diameter());

  // Smooth (diffuse) deposited particle data
  diffusion.smooth( GetVecOfPtrs(m_data), 0, ncomps, m_variables);
}


void
pm_solids::
reset ()
{
  int ncomps = static_cast<int>(m_variables.size());
  int const alpha_p_comp = m_variable_map["alpha"];

  for (int lev(0); lev<get_nlev(); ++lev) {

    for (MFIter mfi(*m_data[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& tbox = mfi.tilebox();

      Array4<Real const> const& alpha = m_data[lev]->const_array(mfi,alpha_p_comp);
      Array4<Real      > const& sdata = m_data[lev]->array(mfi);

      ParallelFor(tbox, ncomps, [alpha, alpha_p_comp, sdata]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { if ((n != alpha_p_comp) && alpha(i,j,k) > 0.)
        { sdata(i,j,k,n) /= alpha(i,j,k); }
      });
    }
  }

  int comp = ncomps;

  m_variable_map["Pp_11"] = comp++;
  m_variable_map["Pp_12"] = comp++;
  m_variable_map["Pp_13"] = comp++;
  m_variable_map["Pp_22"] = comp++;
  m_variable_map["Pp_23"] = comp++;
  m_variable_map["Pp_21"] = m_variable_map["Pp_12"] ;
  m_variable_map["Pp_31"] = m_variable_map["Pp_13"] ;
  m_variable_map["Pp_32"] = m_variable_map["Pp_23"] ;
  m_variable_map["Pp_33"] = comp++;
  m_variable_map["Theta"] = comp++;

  for (int lev(0); lev<get_nlev(); ++lev) {

    MultiFab tmp_MF(m_grids[lev], m_dmap[lev], 1, 0);

    { int const dstcomp = m_variable_map["Pp_11"];
      // Copy Up^2 into m_data
      MultiFab::Copy(*m_data[lev].get(), *m_data[lev].get(), m_variable_map["Up^2"], dstcomp, 1, 0);
      // Compute Up*Up
      MultiFab::Copy(tmp_MF, *m_data[lev].get(), m_variable_map["Up"], 0, 1, 0);
      MultiFab::Multiply(tmp_MF, *m_data[lev].get(), m_variable_map["Up"], 0, 1, 0);
      // Result: Up^2 - Up*Up
      MultiFab::Subtract(*m_data[lev].get(), tmp_MF, 0, dstcomp, 1, 0);
    }

    { int const dstcomp = m_variable_map["Pp_12"];
      // Copy (UpVp) into m_data
      MultiFab::Copy(*m_data[lev].get(), *m_data[lev].get(), m_variable_map["Up*Vp"], dstcomp, 1, 0);
      // Compute Up*Vp
      MultiFab::Copy(tmp_MF, *m_data[lev].get(), m_variable_map["Up"], 0, 1, 0);
      MultiFab::Multiply(tmp_MF, *m_data[lev].get(), m_variable_map["Vp"], 0, 1, 0);
      // Result: (UpVp) - Up*Vp
      MultiFab::Subtract(*m_data[lev].get(), tmp_MF, 0, dstcomp, 1, 0);
    }

    { int const dstcomp = m_variable_map["Pp_13"];
      // Copy (UpWp) into m_data
      MultiFab::Copy(*m_data[lev].get(), *m_data[lev].get(), m_variable_map["Up*Wp"], dstcomp, 1, 0);
      // Compute Up*Wp
      MultiFab::Copy(tmp_MF, *m_data[lev].get(), m_variable_map["Wp"], 0, 1, 0);
      MultiFab::Multiply(tmp_MF, *m_data[lev].get(), m_variable_map["Wp"], 0, 1, 0);
      // Result: (UpWp) - Up*Wp
      MultiFab::Subtract(*m_data[lev].get(), tmp_MF, 0, dstcomp, 1, 0);
    }

    { int const dstcomp = m_variable_map["Pp_22"];
      // Copy Vp^2 into m_data
      MultiFab::Copy(*m_data[lev].get(), *m_data[lev].get(), m_variable_map["Vp^2"], dstcomp, 1, 0);
      // Compute Vp*Vp
      MultiFab::Copy(tmp_MF, *m_data[lev].get(), m_variable_map["Vp"], 0, 1, 0);
      MultiFab::Multiply(tmp_MF, *m_data[lev].get(), m_variable_map["Vp"], 0, 1, 0);
      // Result: Vp^2 - Vp*Vp
      MultiFab::Subtract(*m_data[lev].get(), tmp_MF, 0, dstcomp, 1, 0);
    }

    { int const dstcomp = m_variable_map["Pp_23"];
      // Copy (VpWp) into m_data
      MultiFab::Copy(*m_data[lev].get(), *m_data[lev].get(), m_variable_map["Vp*Wp"], dstcomp, 1, 0);
      // Compute Vp*Wp
      MultiFab::Copy(tmp_MF, *m_data[lev].get(), m_variable_map["Vp"], 0, 1, 0);
      MultiFab::Multiply(tmp_MF, *m_data[lev].get(), m_variable_map["Wp"], 0, 1, 0);
      // Result: (VpWp) - Vp*Wp
      MultiFab::Subtract(*m_data[lev].get(), tmp_MF, 0, dstcomp, 1, 0);
    }

    { int const dstcomp = m_variable_map["Pp_33"];
      // Copy Wp^2 into m_data
      MultiFab::Copy(*m_data[lev].get(), *m_data[lev].get(), m_variable_map["Wp^2"], dstcomp, 1, 0);
      // Compute Wp*Vp
      MultiFab::Copy(tmp_MF, *m_data[lev].get(), m_variable_map["Wp"], 0, 1, 0);
      MultiFab::Multiply(tmp_MF, *m_data[lev].get(), m_variable_map["Wp"], 0, 1, 0);
      // Result: Wp^2 - Wp*Wp
      MultiFab::Subtract(*m_data[lev].get(), tmp_MF, 0, dstcomp, 1, 0);
    }

    { int const dstcomp = m_variable_map["Theta"];
      // Copy Pp_11 into m_data
      MultiFab::Copy(*m_data[lev].get(), *m_data[lev].get(), m_variable_map["Pp_11"], dstcomp, 1, 0);
      // Add Pp_22 and Pp_33
      MultiFab::Add(*m_data[lev].get(), *m_data[lev].get(), m_variable_map["Pp_22"], dstcomp, 1, 0);
      MultiFab::Add(*m_data[lev].get(), *m_data[lev].get(), m_variable_map["Pp_33"], dstcomp, 1, 0);
      // Divide by 3
      m_data[lev]->mult(1.0/3.0, dstcomp, 1, 0);
    }

    m_variables.push_back("Pp_11");
    m_variables.push_back("Pp_12");
    m_variables.push_back("Pp_13");
    m_variables.push_back("Pp_22");
    m_variables.push_back("Pp_23");
    m_variables.push_back("Pp_33");
    m_variables.push_back("Theta");

    m_data[lev]->FillBoundary(m_geom[lev].periodicity());

  }//lev

}
