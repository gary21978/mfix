#ifndef _MFIX_RW_H_
#define _MFIX_RW_H_

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_BoxArray.H>

#ifdef MFIX_MPMD
#include <AMReX_MPMD.H>
#endif

#include <mfix_timer.H>
#include <mfix_pc.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_regions.H>
#include <mfix_monitors.H>

#include <mfix_io_reports_diameter.H>
#include <mfix_mass_balance.H>

class MFIXReadWrite
{
  public:
   MFIXReadWrite (int nlev_in,
                  amrex::Vector<amrex::BoxArray>& grids_in,
                  amrex::Vector<amrex::Geometry>& geom_in,
                  MFIXParticleContainer* pc_in,
                  MFIXFluidPhase& fluid_in,
                  MFIXLevelData& a_leveldata,
                  amrex::Real const& a_therm_p,
                  amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory_in,
                  amrex::Vector<amrex::DistributionMapping>& dmap_in,
                  bool ooo_debug_in,
                  amrex::Vector<std::unique_ptr<amrex::MultiFab>> const& level_sets_in,
                  int levelset_refinement_in,
                  int levelset_pad_in,
                  int levelset_eb_refinement_in,
                  int levelset_eb_pad_in,
                  MFIXSolidsPhase& solids_in,
                  MFIXDEM& dem,
                  MFIXPIC& pic,
                  MFIXReactions& reactions_in,
                  const amrex::Vector<amrex::IntVect>& ref_ratio_in,
                  BCList& bc_list_in,
                  amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& particle_ebfactory_in,
                  std::shared_ptr<MFIXMassBalance> a_mass_balance,
                  MFIXRegions& regions_in,
                  MFIXBoundaryConditions& boundary_conditions_in,
                  const amrex::Vector<amrex::BCRec>& bcrec_velocity_in);

   ~MFIXReadWrite ();

   // Class for specifying a subregion for ParticleContainer::WritePlotFile
   // and add a filter to plot only a selected group of particles
   struct SolidsPlotRegion {
     public:
       SolidsPlotRegion ()
         : m_region_name("")
         , m_region_extents(amrex::RealBox())
         , m_pp_string("")
         , m_plot_names(0)
         , m_h_plot_types(0)
         , m_d_plot_types(0)
         , m_plot_int(-1)
         , m_plot_per_approx(-1.)
         , m_plot_fluid_vars(0)
       {}

       std::string m_region_name;
       amrex::RealBox m_region_extents;

       std::string m_pp_string;

       amrex::Vector<std::string> m_plot_names;
       amrex::Gpu::HostVector<int> m_h_plot_types;
       amrex::Gpu::DeviceVector<int> m_d_plot_types;

       int m_plot_int;
       amrex::Real m_plot_per_approx;
       amrex::Vector<std::string> m_plot_fluid_vars;

       int m_last_solids_plt = -1;

       amrex::Vector<int> m_write_real_comp = amrex::Vector<int>(SoArealData::count,1);
       amrex::Vector<int> m_write_int_comp = amrex::Vector<int>(SoAintData::count,1);
   };

   void set_pc (MFIXParticleContainer* pc_in) {

     this->pc = pc_in;

     // Provide names for all particle data objects. These are used when
     // writing the plot files.
     GetSolidsNames();

     // Set IO flags for particle data for plot files.
     GetSolidsIOPltFlags("mfix", write_real_comp, write_int_comp);

     // Set IO flags for particle sub-region plot files.
     for (int n(0); n < m_solids_plot_regions.size(); ++n) {

       auto& plot_region = m_solids_plot_regions[n];

       GetSolidsIOPltFlags(plot_region.m_pp_string,
           plot_region.m_write_real_comp, plot_region.m_write_int_comp);

     }
   }

   void set_monitors_pc (MFIXParticleContainer* pc_in) {
     m_monitors.set_pc(pc_in);
   }

   void output (int estatus,
                int finish,
                int nstep,
                amrex::Real dt,
                amrex::Real time)
   { pc->output(estatus, finish, nstep, dt, time); }

   void InitIOPltData ();

   void writeNow (MFIXTimer& timer,
                  bool first=false, /* first call*/
                  bool last=false); /* last call*/

   void writeStaticPlotFiles () const;

   void reportGridStats () const;

   void WriteCheckPointFile (std::string& check_file_name,
                             int nstep,
                             amrex::Real dt,
                             amrex::Real time);

   void WritePlotFile (std::string& plot_file_name,
                       int nstep = 0,
                       amrex::Real time = 0.0);

   void WriteSolidsPlotFile (SolidsPlotRegion& plot_region,
                             int nstep = 0,
                             amrex::Real time = 0.0);

   void WriteAscentFile (int nstep,
                         amrex::Real time) const;

#ifdef AMREX_USE_CATALYST
    void RunCatalystAdaptor (int nstep,
                             amrex::Real time) const;

    void init_catalyst_fluid ( amrex::Vector<std::string>& a_names,
                               amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_data ) const;

    void init_catalyst_pc ( amrex::Vector<std::string>& a_names_int,
                            amrex::Vector<std::string>& a_names_real ) const;
#endif

   void WriteParticleAscii (std::string& par_ascii_file_name,
                            int nstep = 0) const;

   void FillPatchVelForVort (int lev, amrex::MultiFab& vel);
   void ComputeVort ();

   void ReportGridStats () const;

   void WriteMassBalanceReport (const amrex::Real new_time);

   void mfix_print_max_vel (int lev,
                            const amrex::Vector<amrex::MultiFab*>& vel_g_in,
                            const amrex::Vector<amrex::MultiFab*>& p_g_in);

   void mfix_print_max_gp (int lev,
                           const amrex::Vector<amrex::MultiFab*>& gp_g_in);

   amrex::IntVect mfix_print_min_epg ();

   void GotoNextLine (std::istream& is);

#ifdef MFIX_MPMD
   void MPMDInit ();
   void MPMDSend (amrex::Real time) const;
   void MPMDFinal () const;
#endif

   int repl_x = 1;
   int repl_y = 1;
   int repl_z = 1;
   bool stop_for_unused_inputs = false;
   bool only_print_grid_report = false;
   std::string restart_file {""};

   int report_mass_balance = 0;
   int mass_balance_report_int = -1;
   int last_mb_report  = -1;
   amrex::Real mass_balance_report_per_approx = -1;
   amrex::Real mass_balance_report_time;

   void setReportTime (amrex::Real const a_time)
   { mass_balance_report_time = a_time; }

#ifdef AMREX_USE_CATALYST
   bool catalyst_on_restart = false;
   std::string catalyst_script {""};
   std::string catalyst_impl {""};
   std::string catalyst_library_path {""};
   bool catalyst_enabled = false;
#endif

   // Variables to simplify checkpoint IO
   // amrex::Vector< amrex::Vector< amrex::MultiFab* > > vectorVars;
   amrex::Vector<std::string> vecVarsName;

   amrex::Vector<amrex::Vector<amrex::MultiFab*>> chkScalarVars;
   amrex::Vector<std::string> chkscaVarsName;

   amrex::Vector<amrex::Vector<amrex::MultiFab*>> chkTVars;
   amrex::Vector<std::string> chkTVarsName;

   amrex::Vector<amrex::Vector<amrex::MultiFab*>> chkSpeciesVars;
   amrex::Vector<std::string> chkSpeciesVarsName;

   const std::string& check_filename () const { return check_file; }

   void Initialize ();

    bool solids_plot_regions () const
    { return m_solids_plot_regions.size() > 0; }

    void set_nlev (int nlev_in) { nlev = nlev_in; }

private:
   // Private members
   int nlev;

   amrex::Vector<amrex::BoxArray>& grids;

   amrex::Vector<amrex::Geometry>& geom;

   MFIXParticleContainer* pc;

   MFIXFluidPhase& fluid;

    MFIXLevelData& m_leveldata;

    MFIXLevelData& leveldata () noexcept { return m_leveldata; }

    MFIXLevelData const& leveldata_const () const noexcept { return m_leveldata; }

    LevelData* leveldata ( int const a_lev ) noexcept {
      AMREX_ASSERT( a_lev < m_leveldata.m_level_data.size() );
      return m_leveldata.m_level_data[a_lev].get();
    }


   amrex::Real const& m_therm_p;

   amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory;
   amrex::Vector<amrex::DistributionMapping>& dmap;

   bool ooo_debug;

   amrex::Vector<std::unique_ptr<amrex::MultiFab>> const& level_sets;

   int levelset_refinement;
   int levelset_pad;
   int levelset_eb_refinement;
   int levelset_eb_pad;

   MFIXSolidsPhase& solids;

   MFIXDEM& m_dem;

   MFIXPIC& m_pic;

   MFIXReactions& reactions;

   const amrex::Vector<amrex::IntVect>& ref_ratio;

   BCList& bc_list;

   amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& particle_ebfactory;

   std::shared_ptr<MFIXMassBalance> m_mass_balance;

   MFIXRegions& regions;

   MFIXBoundaryConditions& m_boundary_conditions;
   const amrex::Vector<amrex::BCRec>& m_bcrec_velocity;

   void readParameters ();

   int test_per_approx(const amrex::Real time,
                       const amrex::Real dt,
                       const amrex::Real per_approx);

   int test_walltime_interval (const MFIXTimer& timer);

   // Flag for writing static particle levelset function in plot file
   bool write_ls = false;
   std::string static_plt_file_ls {"plt_ls"};

   // Flag for writing static EB geometry data in plot file
   bool plt_geom = false;
   std::string static_plt_file_geom {"plt_eb_geom"};

   // minimum digits in filename
   int m_min_digits = 5;

   // Checkpoint file controls
   bool checkpoint_files_output = true;
   std::string check_file {"chk"};

   int check_int = -1;
   amrex::Real check_per_approx = -1.;
   amrex::Real check_walltime_interval = -1.;
   int m_interval_nb = 1;
   //amrex::Real check_per_exact  = -1.;
   int last_chk  = -1;

   // Plot file controls
   std::string plot_file {"plt"};
   bool plotfile_on_restart = false;

   int plot_int = -1;
   amrex::Real plot_per_approx = -1.;
   //amrex::Real plot_per_exact  = -1.;

   // Ascent controls
   bool ascent_on_restart = false;

   int ascent_int = -1;
   amrex::Real ascent_per_approx = -1.;

   // Particle ASCII output files (debug only)
   int par_ascii_int = -1;
   int last_par_ascii  = -1;
   std::string par_ascii_file {"par"};
   amrex::Real par_ascii_per_approx = -1;

   // Last step at which we wrote a plotfile
   int last_plt = -1;

   // Flags for saving fluid data in plot files
   int plt_vel       = 1;
   int plt_epf       = 1;
   int plt_p_g       = 0;
   int plt_ro_g      = 0;
   int plt_MW_g      = 0;
   int plt_trac      = 0;
   int plt_cp_g      = 0;
   int plt_T_g       = 0;
   int plt_h_g       = 0;
   int plt_k_g       = 0;
   int plt_mu_g      = 0;
   int plt_X_gk      = 0;
   int plt_D_gk      = 0;
   int plt_cp_gk     = 0;
   int plt_h_gk      = 0;
   int plt_volfrac   = 0;
   int plt_gradp_g   = 0;
   int plt_vort      = 0;
   int plt_txfr      = 0;
   int plt_chem_txfr = 0;
   int plt_proc      = 0;
   int plt_proc_p    = 0;
   int plt_cost_p    = 0;
   // Geometry data
   int plt_impf      = 0;
   int plt_centroid  = 0;
   int plt_bndryarea = 0;
   int plt_bndrycent = 0;
   int plt_bndrynorm = 0;
   int plt_areafrac  = 0;
   int plt_facecent  = 0;
   int plt_edgecent  = 0;

   // Total number of variables to write in plot file
   int pltVarCount = 0;

   // Flags for saving particle data. By default, we have all flags on,
   // we turn off what we don't want in the init IO routine. This is somewhat
   // different from what we do with the fluid.
   amrex::Vector<int> write_real_comp = amrex::Vector<int>(SoArealData::count,1);
   amrex::Vector<int> write_int_comp = amrex::Vector<int>(SoAintData::count,1);

   amrex::Vector<std::string> real_comp_names;
   amrex::Vector<std::string>  int_comp_names;

   amrex::Vector<SolidsPlotRegion> m_solids_plot_regions;

   // Construct arrays of names for particle data.
   void GetSolidsNames ();

   void GetSolidsIOPltFlags ( std::string const solids_region_prefix,
                              amrex::Vector<int>& a_write_real_comp_out,
                              amrex::Vector<int>& a_write_int_comp_out );

   void WriteCheckHeader (const std::string& name,
                          int nstep,
                          amrex::Real dt,
                          amrex::Real time) const;

   void WriteJobInfo (const std::string& dir) const;

   //! Save variables that don't change to plot file. The idea is that they can
   //! be saved _once_ per simulations and not repeatedly every Nth time step.
   void WriteStaticPlotFileParticleLevelSet (const std::string& plotfilename) const;
   void WriteStaticPlotFileEBGeometry (const std::string& plotfilename) const;

   Monitors m_monitors;

   std::string m_ascent_actions_yaml;

   reports::diameter* m_p_reports_diameter;

#ifdef MFIX_MPMD
   int m_mpmd_int = -1;
   amrex::Real m_mpmd_per_approx = -1.;
   int m_mpmd_this_root = -1;
   int m_mpmd_other_root = -1;

   amrex::Vector<std::unique_ptr<amrex::MPMD::Copier>> m_mpmd_copiers;
   amrex::Vector<std::string> m_mpmd_static_mf_names;
   amrex::Vector<std::string> m_mpmd_mf_names;

   void ConstructMPMDHeader (std::string& header) const;
   void SetMPMDRoots ();
   void SendMPMDIntFlags (bool end_flag) const;
   void SendMPMDReals (amrex::Real time) const;
   void SendMPMDMultiFabs () const;
#endif
};

#endif
