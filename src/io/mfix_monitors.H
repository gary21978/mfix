#ifndef MFIX_MONITORS_H_
#define MFIX_MONITORS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBFArrayBox.H>

#include <mfix_pc.H>
#include <mfix_regions.H>
#include <mfix_leveldata.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_run_on.H>

#include <iostream>
#include <filesystem>
#include <map>
#include <string>
#include <array>
#include <fstream>


// forward declaration
class Monitor;


// Container for all the possible monitors
class Monitors {
  public:
    // Constructor
    explicit Monitors () = default;

    // Destructor
    ~Monitors ();

    void
    initialize (MFIXRegions const& regions,
                MFIXLevelData& a_leveldata,
                amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                MFIXFluidPhase& fluid,
                MFIXParticleContainer* pc,
                amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& particle_ebfactory,
                MFIXSolidsPhase& solids);

    Monitor& get (const int i)
    {
      AMREX_ASSERT(i >=0 && i < m_monitors.size());
      return *m_monitors[i];
    }

    int size () const
    { return m_monitors.size(); }

    void set_pc (MFIXParticleContainer* pc);

  private:
    amrex::Vector<Monitor*> m_monitors;
};


// Single monitor class
class Monitor {

  public:
    // Constructor
    Monitor (amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
             amrex::RealBox const& region,
             amrex::RealBox const* plane = nullptr);

    // Constructor
    Monitor (const std::array<std::string,2> specs,
             const std::string pparse,
             amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
             MFIXRegions const& regions);

    // Destructor
    virtual ~Monitor () = default;

    const int & plot_int () const
    { return m_plot_int; }

    const amrex::Real& plot_per_approx () const
    { return m_plot_per_approx; }

    void write_csv (const amrex::Real& time,
                    const amrex::Real& dt);

    virtual
    void set_pc (MFIXParticleContainer* pc) = 0;

  protected:
    // Setup the monitor
    virtual
    void initialize ();

    // Initialization flag
    int is_initialized;

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) = 0;

    // The function to count the number of cells or particles
    virtual
    void count () = 0;

    // Number of levels
    int m_nlev;

    // Set variables, components and indexes
    virtual
    void setup_variables () = 0;

    // Monitor specifications
    const std::array<std::string,2> m_specs;

    // Data structure where to store the results of the monitoring action
    amrex::Vector<amrex::Vector<amrex::Real>> m_monitoring_results;

    // Count cells or particles if the flag is ON
    int m_count_flag;
    amrex::Vector<amrex::Real> m_counts;

    // Reference to the EBFactory
    amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& m_ebfactory;

    // Reference to the mfix regions object
    MFIXRegions const* m_regions;

    // Filename for the output file
    std::filesystem::path m_inputs_full_path;
    std::vector<std::filesystem::path> m_full_paths;
    std::string m_filename;

    // Timestep number for outputting
    int m_plot_int;

    // Simulation time for outputting
    amrex::Real m_plot_per_approx;

    // Monitor output file openmode
    std::string m_openmode;

    // Output formatting options
    int m_setw;
    std::string m_setfill;
    int m_setprecision;
    std::string m_formatflag;

    // Apply formatting options
    std::ostream& setformat(std::ostream& output) const;

    // The name of the region for the monitor
    std::string m_region_name;

    // The RealBox corresponding to the monitored region
    const amrex::RealBox* m_region;

    // The variable to plot
    amrex::Vector<std::string> m_variables;
    amrex::Vector<std::string> m_input_variables;

    // The Box corresponding to the monitored region
    amrex::Vector<amrex::Box> m_boxes;

    // Area normal vector direction: {0: x_direction, 1: y_direction, 2: z_direction}
    int m_direction;

    // The name of the region for the flow rate plane
    std::string m_plane_name;

    // The RealBox needed for Lagrangian monitors of type FlowRate
    const amrex::RealBox* m_plane;

    // Function for computing the Box out of the RealBox
    amrex::Box calc_box (const amrex::Geometry&, const amrex::RealBox&) const;

    // Function for computing the Box out of the RealBox
    amrex::RealBox calc_realbox (const amrex::Geometry&, const amrex::Box&) const;

    // Function for computing the Box out of the RealBox
    amrex::RealBox realboxes_intersection (const amrex::RealBox&, const amrex::RealBox&) const;
};


namespace EulerianMonitor {

class BaseMonitor: public Monitor {

  public:
    // Constructor
    BaseMonitor ( MFIXLevelData& a_leveldata,
                  amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                  amrex::RealBox const& region,
                  amrex::RealBox const* plane = nullptr)
      : Monitor(ebfactory, region, plane)
      , m_leveldata(a_leveldata)
    {
      AMREX_ALWAYS_ASSERT(m_nlev > 0);

      m_epsilon.resize(m_nlev, nullptr);
      m_density.resize(m_nlev, nullptr);
      m_velocity.resize(m_nlev, nullptr);

      for (int lev(0); lev < m_nlev; ++lev) {
        m_epsilon[lev] = leveldata().epf(lev);
        m_density[lev] = leveldata().rho(lev);
        m_velocity[lev] = leveldata().vel(lev);
      }
    }

    // Constructor
    BaseMonitor (const std::array<std::string,2> specs,
                 const std::string& pparse,
                 MFIXLevelData& a_leveldata,
                 amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                 MFIXFluidPhase& fluid,
                 MFIXRegions const& regions)
      : Monitor(specs, pparse, ebfactory, regions)
      , m_leveldata(a_leveldata)
      , m_fluid(&fluid)
    {}

    virtual
    void set_pc (MFIXParticleContainer*) override
    {}

    // Template class for generic computation
    template <typename D, typename Op, typename F>
    typename D::Type
    apply (const int lev,
           const amrex::MultiFab* mf,
           const int component,
           D& reduce_data,
           Op& reduce_op,
           F&& f,
           typename D::Type& default_values)
    {
      using ReduceTuple = typename D::Type;

      ReduceTuple result(default_values);

      const amrex::MultiFab* mf_copy = mf;

      int mf_converted = convert_mf_if_needed(mf, lev);
      int epg_converted = convert_mf_if_needed(m_epsilon[lev], lev);
      int rog_converted = convert_mf_if_needed(m_density[lev], lev);
      int velg_converted = convert_mf_if_needed(m_velocity[lev], lev);

      int vel_Comp(0);

      if (m_velocity[lev]) {
        if (m_velocity[lev]->nComp() == 3) {
          vel_Comp = m_direction;
        } else if (m_velocity[lev]->nComp() == 1) {
          vel_Comp = 0;
        } else {
          amrex::Abort("Error: velocity MultiFab has wrong number of components");
        }
      }

      amrex::BoxArray box_array = this->get_box_array(lev);
      const amrex::DistributionMapping& d_map = this->get_d_map(lev);

      const amrex::FabArray<amrex::EBCellFlagFab>& flags = m_ebfactory[lev]->getMultiEBCellFlagFab();

      bool was_there_intersection(false);

      for (amrex::MFIter mfi(box_array, d_map, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = mfi.tilebox();

        const amrex::EBCellFlagFab& flags_fab = flags[mfi];

        bool single_valued_fab(false);
        if (flags_fab.getType() == amrex::FabType::singlevalued) {
          single_valued_fab = true;
        }

        if (bx.intersects(m_boxes[lev])) {

          const amrex::Box l_box = bx & m_boxes[lev];

          amrex::FArrayBox farray_box(l_box);
          if (!single_valued_fab)
            farray_box.setVal<run_on>(1.);

          AMREX_ASSERT(l_box.ok());

          const int grow_nb = this->get_box_grow_nb();

          // Only if fab is not covered
          if (flags_fab.getType(amrex::grow(l_box, grow_nb)) != amrex::FabType::covered) {
            was_there_intersection = true;

            amrex::Array4<amrex::Real const> dummy_arr;

            auto const& mf_arr = mf ? mf->const_array(mfi,component) : dummy_arr;
            auto const& fraction_arr = single_valued_fab ? get_fraction_arr(mfi,lev) : farray_box.const_array();
            auto const& epsilon_arr = m_epsilon[lev] ? m_epsilon[lev]->const_array(mfi) : dummy_arr;
            auto const& density_arr = m_density[lev] ? m_density[lev]->const_array(mfi) : dummy_arr;
            auto const& velocity_arr = m_velocity[lev] ? m_velocity[lev]->const_array(mfi,vel_Comp) : dummy_arr;
            auto const& flags_arr = flags_fab.const_array();

            isCovered is_covered = m_is_covered;

            reduce_op.eval(l_box, reduce_data, [mf_arr,fraction_arr,epsilon_arr,
                density_arr,velocity_arr,f,flags_arr,default_values,is_covered]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
              if (!is_covered(flags_arr, amrex::IntVect(i,j,k))) {

                return f(mf_arr, fraction_arr, epsilon_arr, density_arr,
                         velocity_arr, amrex::IntVect(i,j,k));
              } else {

                return default_values;
              }
            });
          } // if l_box is ok
        } // if fab not covered
      } // MFIter loop

      if (mf && mf_converted) { delete mf; mf = mf_copy; }
      if (m_epsilon[lev] && epg_converted) { delete m_epsilon[lev]; m_epsilon[lev] = leveldata().epf(lev); }
      if (m_density[lev] && rog_converted) { delete m_density[lev]; m_density[lev] = leveldata().rho(lev); }
      if (m_velocity[lev] && velg_converted) { delete m_velocity[lev]; m_velocity[lev] = leveldata().vel(lev); }

      result = was_there_intersection ? reduce_data.value(reduce_op) : default_values;

      return result;
    }

    // Description: Class that handles different situations when it is needed to
    // determine if a cell is covered or not
    class isCovered {

      public:
        AMREX_GPU_HOST_DEVICE
        isCovered ()
          : m_ix_type(-1)
          , m_direction(-1)
        {}

        AMREX_GPU_HOST_DEVICE
        isCovered (const int ix_type,
                   const int direction = -1)
          : m_ix_type(ix_type)
          , m_direction(direction)
        {}

        AMREX_GPU_HOST_DEVICE
        void set_direction (const int direction) { m_direction = direction; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int operator() (const amrex::Array4<amrex::EBCellFlag const>& flags_arr,
                        const amrex::IntVect& ijk) const
        {
          if (m_ix_type == 0) {

            return flags_arr(ijk).isCovered();

          } else if (m_ix_type == 1) {

            amrex::IntVect ijk_minus(ijk);
            ijk_minus[m_direction] -= 1;
            return flags_arr(ijk).isCovered() && flags_arr(ijk_minus).isCovered();

          } else {

            amrex::Abort("Wrong ixType");
          }

          return 0;
        }

      private:
        int m_ix_type;
        int m_direction;
    };

    // The function to count the number of cells
    virtual
    void count () override;

  protected:
    // Set variables, components and indexes
    virtual
    void setup_variables () override;

    // The variables
    amrex::Vector<amrex::Vector<const amrex::MultiFab*>> m_mf;
    amrex::Vector<amrex::Vector<const amrex::Real*>> m_scalars;

    // The components
    amrex::Vector<int> m_components;

    // Convert MultiFabs when they're not of the correct type
    virtual
    int convert_mf_if_needed (const amrex::MultiFab*& mf,
                              const int lev) = 0;

    // Setup the monitor
    virtual
    void initialize () override
    { Monitor::initialize(); }

    MFIXLevelData& m_leveldata;

    MFIXLevelData& leveldata () noexcept { return m_leveldata; }

    LevelData* leveldata ( int const a_lev ) noexcept {
      AMREX_ASSERT( a_lev < m_leveldata.m_level_data.size() );
      return m_leveldata.m_level_data[a_lev].get();
    }

    MFIXLevelData const& leveldata_const () const noexcept { return m_leveldata; }

    amrex::Vector<const amrex::MultiFab*> m_epsilon;
    amrex::Vector<const amrex::MultiFab*> m_density;
    amrex::Vector<const amrex::MultiFab*> m_velocity;

    amrex::Vector<amrex::MultiFab> m_ones;

    // Reference to the FluidPhase class
    MFIXFluidPhase* m_fluid;

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override = 0;

    // Get the number of cells to grow the tile box when checking for covered FAB
    virtual
    int get_box_grow_nb () const = 0;

    // Check the box is good
    virtual
    void check_boxes_are_ok () const = 0;

    // Check a given MultiFab is ok
    virtual
    bool check_mf_is_ok (const amrex::MultiFab& mf) const = 0;

    // A dummy Array4 needed as a return argument for get_fraction_arr method
    amrex::Array4<amrex::Real const> m_dummy_arr;

    // Return the areafrac or volfrac Array4 object
    virtual
    amrex::Array4<amrex::Real const> get_fraction_arr (amrex::MFIter&,
                                                       const int lev) const = 0;

    // Return the BoxArray for building the MFIter
    virtual
    amrex::BoxArray get_box_array (const int lev) const = 0;

    // Return the DistributionMapping for building the MFIter
    virtual
    const amrex::DistributionMapping& get_d_map (const int lev) const = 0;

    isCovered m_is_covered;
};


class UniformScalarField final: public BaseMonitor {

  public:
    // Constructor
    UniformScalarField ( MFIXLevelData& a_leveldata,
                         amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                         amrex::RealBox const& region)
      : BaseMonitor(a_leveldata, ebfactory, region)
    {}

    // Constructor
    UniformScalarField (const std::array<std::string,2> specs,
                        const std::string& pparse,
                        MFIXLevelData& a_leveldata,
                        amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                        MFIXFluidPhase& fluid,
                        MFIXRegions const& regions)
      : BaseMonitor(specs, pparse, a_leveldata, ebfactory, fluid, regions)
    {}

    // Returns the value of the field quantity in the selected region
    amrex::Vector<amrex::Real> value (const amrex::Vector<const amrex::Real*>& scalars);

  protected:
    virtual
    int convert_mf_if_needed (const amrex::MultiFab*&,
                              const int) override { return 0; }

    // Setup the monitor
    virtual
    void initialize () override { BaseMonitor::initialize(); }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override;

    // Get the number of cells to grow the tile box when checking for covered FAB
    virtual
    int get_box_grow_nb () const override { return 0; }

    // Check the box is good
    virtual
    void check_boxes_are_ok () const override {}

    // Check a given MultiFab is ok
    virtual
    bool check_mf_is_ok (const amrex::MultiFab&) const override { return true; }

    // Return the volfrac Array4 object
    virtual
    amrex::Array4<amrex::Real const> get_fraction_arr (amrex::MFIter& mfi,
                                                       const int lev) const override
    { return m_ebfactory[lev]->getVolFrac().const_array(mfi); }

    // Return the BoxArray for building the MFIter
    virtual
    amrex::BoxArray get_box_array (const int lev) const override
    { return m_ebfactory[lev]->boxArray(); }

    // Return the DistributionMapping for building the MFIter
    virtual
    const amrex::DistributionMapping& get_d_map (const int lev) const override
    { return m_ebfactory[lev]->DistributionMap(); }
};


class PointRegion final: public BaseMonitor {

  public:
    // Constructor
    PointRegion ( MFIXLevelData& a_leveldata,
                  amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                  amrex::RealBox const& region)
      : BaseMonitor(a_leveldata, ebfactory, region)
    {
      m_is_covered = isCovered(0);

      this->check_boxes_are_ok();
      m_direction = 0;

      this->set_point_coords();
    }

    // Constructor
    PointRegion (const std::array<std::string,2> specs,
                 const std::string& pparse,
                 MFIXLevelData& a_leveldata,
                 amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                 MFIXFluidPhase& fluid,
                 MFIXRegions const& regions)
      : BaseMonitor(specs, pparse, a_leveldata, ebfactory, fluid, regions)
    {
      m_is_covered = isCovered(0);
    }

    // Returns the value of the field quantity in the selected region
    amrex::Vector<amrex::Real> value (const int lev,
                                      const amrex::Vector<const amrex::MultiFab*>& mf,
                                      const amrex::Vector<int>& components);

  protected:
    virtual
    int convert_mf_if_needed (const amrex::MultiFab*& mf,
                              const int lev) override;

    // Setup the monitor
    virtual
    void initialize () override
    {
      BaseMonitor::initialize();

      this->check_boxes_are_ok();
      m_direction = 0;

      this->set_point_coords();
    }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override;

    // Set the m_point coordinates from the m_region RealBox read from inputs
    void set_point_coords ();

    // The point at which the MultiFabs will be evaluated
    amrex::RealVect m_point;

    // Get the number of cells to grow the tile box when checking for covered FAB
    virtual
    int get_box_grow_nb () const override
    { return 0; }

    // Check the box is good
    virtual
    void check_boxes_are_ok () const override;

    // Check a given MultiFab is ok
    virtual
    bool check_mf_is_ok (const amrex::MultiFab& mf) const override;

    // Return the volfrac Array4 object
    virtual
    amrex::Array4<amrex::Real const> get_fraction_arr (amrex::MFIter& mfi,
                                                       const int lev) const override
    {
      return m_ebfactory[lev]->getVolFrac().const_array(mfi);
    }

    // Return the BoxArray for building the MFIter
    virtual
    amrex::BoxArray get_box_array (const int lev) const override
    { return m_ebfactory[lev]->boxArray(); }

    // Return the DistributionMapping for building the MFIter
    virtual
    const amrex::DistributionMapping& get_d_map (const int lev) const override
    { return m_ebfactory[lev]->DistributionMap(); }

  private:
    // Template class for generic computation
    template <typename D, typename Op, typename F>
    void
    apply (const int, const amrex::MultiFab*, const int, D&, Op&, F&&, typename D::Type&)
    {}
};


class AreaMonitor: public BaseMonitor {

  public:
    // Constructor
    AreaMonitor ( MFIXLevelData& a_leveldata,
                  amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                  amrex::RealBox const& region,
                  amrex::RealBox const* plane = nullptr)
      : BaseMonitor(a_leveldata, ebfactory, region, plane)
    {
      m_is_covered = isCovered(1);

      this->check_boxes_are_ok();
      this->set_direction();

      m_is_covered.set_direction(m_direction);

      m_epsilon.clear();
      m_epsilon.resize(m_nlev, nullptr);

      m_density.clear();
      m_density.resize(m_nlev, nullptr);

      m_velocity.clear();
      m_velocity.resize(m_nlev, nullptr);

      // Convert m_epsilon, m_density, and m_velocity if needed
      for (int lev(0); lev < m_nlev; ++lev) {
        m_epsilon[lev] = leveldata().epf(lev);
        m_density[lev] = leveldata().rho(lev);
        m_velocity[lev] = leveldata().vel(lev);
      }
    }

    // Constructor
    AreaMonitor (const std::array<std::string,2> specs,
                 const std::string& pparse,
                 MFIXLevelData& a_leveldata,
                 amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                 MFIXFluidPhase& fluid,
                 MFIXRegions const& regions)
      : BaseMonitor(specs, pparse, a_leveldata, ebfactory, fluid, regions)
    {
      m_is_covered = isCovered(1);
    }

  protected:
    virtual
    int convert_mf_if_needed (const amrex::MultiFab*& mf,
                              const int lev) override;

    // Setup the monitor
    virtual
    void initialize () override
    {
      BaseMonitor::initialize();

      this->check_boxes_are_ok();
      this->set_direction();

      m_is_covered.set_direction(m_direction);
    }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override = 0;

    // Get the number of cells to grow the tile box when checking for covered FAB
    virtual
    int get_box_grow_nb () const override
    { return 1; }

    // Check the box is good
    virtual
    void check_boxes_are_ok () const override;

    // Check a given MultiFab is face centered along m_direction
    virtual
    bool check_mf_is_ok (const amrex::MultiFab& mf) const override;

    // Function to set the area normal vector direction and check the box is 2D
    virtual void set_direction ();

    // Return the areafrac Array4 object
    virtual
    amrex::Array4<amrex::Real const> get_fraction_arr (amrex::MFIter& mfi,
                                                       const int lev) const override
    {
      return m_ebfactory[lev]->getAreaFrac()[m_direction]->const_array(mfi);
    }

    // Return the BoxArray for building the MFIter
    virtual
    amrex::BoxArray get_box_array (const int lev) const override
    {
      amrex::BoxArray edge_ba = m_ebfactory[lev]->boxArray();
      edge_ba.surroundingNodes(m_direction);
      return edge_ba;
    }

    // Return the DistributionMapping for building the MFIter
    virtual
    const amrex::DistributionMapping& get_d_map (const int lev) const override
    { return m_ebfactory[lev]->DistributionMap(); }
};


class VolumeMonitor: public BaseMonitor {

  public:
    // Constructor
    VolumeMonitor ( MFIXLevelData& a_leveldata,
                    amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                    amrex::RealBox const& region)
      : BaseMonitor(a_leveldata, ebfactory, region)
    {
      m_is_covered = isCovered(0);

      this->check_boxes_are_ok();
      m_direction = 0;

      m_epsilon.clear();
      m_epsilon.resize(m_nlev, nullptr);

      m_density.clear();
      m_density.resize(m_nlev, nullptr);

      m_velocity.clear();
      m_velocity.resize(m_nlev, nullptr);

      // Convert m_epsilon, m_density, and m_velocity if needed
      for (int lev(0); lev < m_nlev; ++lev) {
        m_epsilon[lev] = leveldata().epf(lev);
        m_density[lev] = leveldata().rho(lev);
        m_velocity[lev] = leveldata().vel(lev);
      }
    }

    // Constructor
    VolumeMonitor (const std::array<std::string,2> specs,
                   const std::string& pparse,
                   MFIXLevelData& a_leveldata,
                   amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                   MFIXFluidPhase& fluid,
                   MFIXRegions const& regions)
      : BaseMonitor(specs, pparse, a_leveldata, ebfactory, fluid, regions)
    {
      m_is_covered = isCovered(0);
    }

  protected:
    virtual
    int convert_mf_if_needed (const amrex::MultiFab*& mf,
                              const int lev) override;

    // Setup the monitor
    virtual
    void initialize () override
    {
      BaseMonitor::initialize();

      this->check_boxes_are_ok();
      m_direction = 0;
    }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override = 0;

    // Get the number of cells to grow the tile box when checking for covered FAB
    virtual
    int get_box_grow_nb () const override
    { return 0; }

    // Check the box is good
    virtual
    void check_boxes_are_ok () const override;

    // Check a given MultiFab is cell centered
    virtual
    bool check_mf_is_ok (const amrex::MultiFab& mf) const override;

    // Return the volfrac Array4 object
    virtual
    amrex::Array4<amrex::Real const> get_fraction_arr (amrex::MFIter& mfi,
                                                       const int lev) const override
    { return m_ebfactory[lev]->getVolFrac().const_array(mfi); }

    // Return the BoxArray for building the MFIter
    virtual
    amrex::BoxArray get_box_array (const int lev) const override
    { return m_ebfactory[lev]->boxArray(); }

    // Return the DistributionMapping for building the MFIter
    virtual
    const amrex::DistributionMapping& get_d_map (const int lev) const override
    { return m_ebfactory[lev]->DistributionMap(); }
};


class AreaRegion final: public AreaMonitor {

  public:
    // Constructor
    AreaRegion ( MFIXLevelData& a_leveldata,
                 amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                 amrex::RealBox const& region)
      : AreaMonitor(a_leveldata, ebfactory, region)
    {
      m_is_covered = isCovered(1);

      m_is_covered.set_direction(m_direction);
    }

    // Constructor
    AreaRegion (const std::array<std::string,2> specs,
                const std::string& pparse,
                MFIXLevelData& a_leveldata,
                amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                MFIXFluidPhase& fluid,
                MFIXRegions const& regions)
      : AreaMonitor(specs, pparse, a_leveldata, ebfactory, fluid, regions)
    {
      m_is_covered = isCovered(1);
    }

    // The sum is computed by summing all the values of the field quantity in
    // the selected region
    amrex::Vector<amrex::Real> sum (const int lev,
                                    const amrex::Vector<const amrex::MultiFab*>& mf,
                                    const amrex::Vector<int>& components);

    // Minimum value of the field quantity in the selected region
    amrex::Vector<amrex::Real> min (const int lev,
                                    const amrex::Vector<const amrex::MultiFab*>& mf,
                                    const amrex::Vector<int>& components);

    // Maximum value of the field quantity in the selected region
    amrex::Vector<amrex::Real> max (const int lev,
                                    const amrex::Vector<const amrex::MultiFab*>& mf,
                                    const amrex::Vector<int>& components);

    // Average value of the field quantity in the selected region
    amrex::Vector<amrex::Real> average (const int lev,
                                        const amrex::Vector<const amrex::MultiFab*>& mf,
                                        const amrex::Vector<int>& components);

    // Standard deviation of the field quantity in the selected region
    amrex::Vector<amrex::Real> stddev (const int lev,
                                       const amrex::Vector<const amrex::MultiFab*>& mf,
                                       const amrex::Vector<int>& components);

  protected:
    // Setup the monitor
    virtual
    void initialize () override
    {
      AreaMonitor::initialize();

      m_is_covered.set_direction(m_direction);
    }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override;
};


class VolumeRegion final: public VolumeMonitor {

  public:
    // Constructor
    VolumeRegion ( MFIXLevelData& a_leveldata,
                   amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                   amrex::RealBox const& region)
      : VolumeMonitor(a_leveldata, ebfactory, region)
    {
      m_is_covered = isCovered(0);
    }

    // Constructor
    VolumeRegion (const std::array<std::string,2> specs,
                  const std::string& pparse,
                  MFIXLevelData& a_leveldata,
                  amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                  MFIXFluidPhase& fluid,
                  MFIXRegions const& regions)
      : VolumeMonitor(specs, pparse, a_leveldata, ebfactory, fluid, regions)
    {
      m_is_covered = isCovered(0);
    }

    // The sum is computed by summing all the values of the field quantity in
    // the selected region
    amrex::Vector<amrex::Real> sum (const int lev,
                                    const amrex::Vector<const amrex::MultiFab*>& mf,
                                    const amrex::Vector<int>& components);

    // Minimum value of the field quantity in the selected region
    amrex::Vector<amrex::Real> min (const int lev,
                                    const amrex::Vector<const amrex::MultiFab*>& mf,
                                    const amrex::Vector<int>& components);

    // Maximum value of the field quantity in the selected region
    amrex::Vector<amrex::Real> max (const int lev,
                                    const amrex::Vector<const amrex::MultiFab*>& mf,
                                    const amrex::Vector<int>& components);

    // Average value of the field quantity in the selected region
    amrex::Vector<amrex::Real> average (const int lev,
                                        const amrex::Vector<const amrex::MultiFab*>& mf,
                                        const amrex::Vector<int>& components);

    // Standard deviation of the field quantity in the selected region
    amrex::Vector<amrex::Real> stddev (const int lev,
                                       const amrex::Vector<const amrex::MultiFab*>& mf,
                                       const amrex::Vector<int>& components);

  protected:
    // Setup the monitor
    virtual
    void initialize () override
    { VolumeMonitor::initialize(); }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override;
};


class SurfaceIntegral final: public AreaMonitor {

  public:
    // Constructor
    SurfaceIntegral ( MFIXLevelData& a_leveldata,
                      amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                      amrex::RealBox const& region,
                      amrex::RealBox const* plane)
      : AreaMonitor(a_leveldata, ebfactory, region, plane)
    {
      m_is_covered = isCovered(1);
      set_dA();
      m_is_covered.set_direction(m_direction);
    }

    // Constructor
    SurfaceIntegral (const std::array<std::string,2> specs,
                     const std::string& pparse,
                     MFIXLevelData& a_leveldata,
                     amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                     MFIXFluidPhase& fluid,
                     MFIXRegions const& regions)
      : AreaMonitor(specs, pparse, a_leveldata, ebfactory, fluid, regions)
    {
      m_is_covered = isCovered(1);
    }

    // Area of the selected region is computed by summing the areas of the
    // facets that define the surface
    amrex::Vector<amrex::Real> area (const int lev);

    // The area-weighted average is computed by dividing the summation of the
    // product of the selected variable and facet area by the total area of the
    // region
    amrex::Vector<amrex::Real> area_weighted_average (const int lev,
                                                      const amrex::Vector<const amrex::MultiFab*>& mf,
                                                      const amrex::Vector<int>& components);

    // The flow rate of a field variable through a surface is computed by
    // summing the product of the phase volume fraction, desity, the selected
    // field variable, phase velocity normal to the facet v_n, and the facet
    // area
    amrex::Vector<amrex::Real> flow_rate (const int lev,
                                          const amrex::Vector<const amrex::MultiFab*>& mf,
                                          const amrex::Vector<int>& components);

    // The mass flow rate through a surface is computed by summing the product
    // of the phase volume fraction, density, phase velocity normal to the facet
    // v_n, and the facet area
    amrex::Vector<amrex::Real> mass_flow_rate (const int lev);

    // The mass flow rate through a surface is computed by summing the
    // product of the phase volume fraction, density, velocity normal to the
    // facet v_n, and the facet area
    amrex::Vector<amrex::Real> mass_weighted_average (const int lev,
                                                      const amrex::Vector<const amrex::MultiFab*>& mf,
                                                      const amrex::Vector<int>& components);

    // The volume flow rate through a surface is computed by summing the product
    // of the phase volume fraction, phase velocity normal to the facet v_n, and
    // the facet area
    amrex::Vector<amrex::Real> volume_flow_rate (const int lev);

  protected:
    // Setup the monitor
    virtual
    void initialize () override
    {
      AreaMonitor::initialize();

      this->set_dA();

      m_is_covered.set_direction(m_direction);
    }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override;

    // Function to set facets area
    void set_dA ();

    // Facets area
    amrex::Vector<amrex::Real> m_dA;
};


class VolumeIntegral final: public VolumeMonitor {

  public:
    // Constructor
    VolumeIntegral ( MFIXLevelData& a_leveldata,
                     amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                     amrex::RealBox const& region)
      : VolumeMonitor(a_leveldata, ebfactory, region)
    {
      m_is_covered = isCovered(0);
      set_dV();
    }

    // Constructor
    VolumeIntegral (const std::array<std::string,2> specs,
                    const std::string& pparse,
                    MFIXLevelData& a_leveldata,
                    amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                    MFIXFluidPhase& fluid,
                    MFIXRegions const& regions)
      : VolumeMonitor(specs, pparse, a_leveldata, ebfactory, fluid, regions)
    {
      m_is_covered = isCovered(0);
    }

    // The volume is computed by summing all of the cell volumes in the selected
    // region
    amrex::Vector<amrex::Real> volume (const int lev);

    // The volume integral is computed by summing the product of the selected
    // field variable and the cell volume
    amrex::Vector<amrex::Real> volume_integral (const int lev,
                                                const amrex::Vector<const amrex::MultiFab*>& mf,
                                                const amrex::Vector<int>& components);

    // The volume integral is computed by summing the product of the selected
    // field variable and the cell volume
    amrex::Vector<amrex::Real> volume_weighted_integral (const int lev,
                                                         const amrex::Vector<const amrex::MultiFab*>& mf,
                                                         const amrex::Vector<int>& components);

    // The volume-weighted average is computed by summing the product of phase
    // volume fraction, density, selected field variable, and cell volume
    amrex::Vector<amrex::Real> volume_weighted_average (const int lev,
                                                        const amrex::Vector<const amrex::MultiFab*>& mf,
                                                        const amrex::Vector<int>& components);

    // The mass-weighted integral is computed by summing the product of phase
    // volume fraction, density, selected field variable, and cell volume
    amrex::Vector<amrex::Real> mass_weighted_integral (const int lev,
                                                       const amrex::Vector<const amrex::MultiFab*>& mf,
                                                       const amrex::Vector<int>& components);

    // The mass-weighted average is computed by dividing the sum of the product
    // of phase volume fraction, density, selected field variable, and cell
    // volume by the summation of the product of the phase volume fraction,
    // density, and cell volume
    amrex::Vector<amrex::Real> mass_weighted_average (const int lev,
                                                      const amrex::Vector<const amrex::MultiFab*>& mf,
                                                      const amrex::Vector<int>& components);

  protected:
    // Setup the monitor
    virtual
    void initialize () override
    {
      VolumeMonitor::initialize();

      this->set_dV();
    }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override;

    // Function to set cells volume
    void set_dV ();

    // Cells volume
    amrex::Vector<amrex::Real> m_dV;
};

} // end namespace EulerianMonitor


namespace LagrangianMonitor {

// General property class
class BaseMonitor: public Monitor {

  public:
    using ParticleTileData = amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>;

    // Constructor
    BaseMonitor (MFIXParticleContainer* pc,
                 amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                 amrex::RealBox const& region,
                 amrex::RealBox const* plane = nullptr)
      : Monitor(ebfactory, region, plane)
      , m_pc(pc)
    {
      this->check_realbox_is_ok();
      m_direction = 0;
    }

    // Constructor
    BaseMonitor (const std::array<std::string,2> specs,
                 const std::string& pparse,
                 MFIXParticleContainer* pc,
                 amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                 MFIXSolidsPhase& solids,
                 MFIXRegions const& regions)
      : Monitor(specs, pparse, ebfactory, regions)
      , m_pc(pc)
      , m_solids(&solids)
    {}

    virtual
    void set_pc (MFIXParticleContainer* pc) override
    { m_pc = pc; }

    // Template class for generic computation
    template <typename D, typename Op, typename F>
    typename D::Type
    apply (const int lev,
           const int index,
           D& reduce_data,
           Op& reduce_op,
           F&& f,
           typename D::Type default_values)
    {
      using ReduceTuple = typename D::Type;

      ReduceTuple result(default_values);

      bool was_there_intersection(false);

      for (MFIXParticleContainer::MFIXParIter pti(*m_pc, lev); pti.isValid(); ++pti) {

        amrex::Box bx = pti.tilebox();
        amrex::RealBox tile_region = calc_realbox(m_ebfactory[lev]->Geom(), bx);

        if (tile_region.intersects(*m_region)) {

          const amrex::RealBox l_region = realboxes_intersection(tile_region, *m_region);

          AMREX_ASSERT(l_region.ok());
          was_there_intersection = true;

          const int np = pti.numRealParticles();

          auto& ptile = pti.GetParticleTile();
          const ParticleTileData& ptile_data = ptile.getParticleTileData();

          reduce_op.eval(np, reduce_data, [ptile_data,f,l_region,default_values,index]
            AMREX_GPU_DEVICE (int i) -> ReduceTuple
          {
            const amrex::RealVect& pos = ptile_data.getSuperParticle(i).pos();

            if (l_region.contains(pos)) {

              return f(ptile_data, index, i);

            } else {

              return default_values;
            }
          });
        }
      }

      result = was_there_intersection ? reduce_data.value() : default_values;

      return result;
    }

    class GetParticleValue {
      public:
        static constexpr int ones = -1;
        static constexpr int volume = -2;
        static constexpr int mass = -3;
        static constexpr int k_energy = -4;

        AMREX_GPU_HOST_DEVICE
        GetParticleValue (const runtimeRealData& rrData)
          : m_rrData_count(rrData.count)
        {}

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::Real operator() (const ParticleTileData& ptile_data,
                                const int index,
                                const int i) const
        {
          const auto& particle = ptile_data.getSuperParticle(i);

          if (index == ones) {

            return 1.;

          } else if (index == volume) {

            return SoArealData::volume(particle.rdata(SoArealData::radius));

          } else if (index == mass) {

            const amrex::Real p_radius = particle.rdata(SoArealData::radius);
            return particle.rdata(SoArealData::density) * SoArealData::volume(p_radius);

          } else if (index == k_energy) {

            const amrex::Real p_radius = particle.rdata(SoArealData::radius);
            const amrex::Real p_density = particle.rdata(SoArealData::density);
            const amrex::Real p_mass = p_density * SoArealData::volume(p_radius);
            const amrex::RealVect p_velocity(AMREX_D_DECL(particle.rdata(SoArealData::velx),
                                                          particle.rdata(SoArealData::vely),
                                                          particle.rdata(SoArealData::velz)));

            return .5 * p_mass * p_velocity.radSquared();

          } else if ((index >= 0) && (index <= 2)) {

            return particle.pos(index);

          } else if (index == 3) {

            return particle.id();

          } else if (index == 4) {

            return particle.cpu();

          } else if (index < (5+SoArealData::count)) {

            const int idx = index-5;
            return particle.rdata(idx);

          } else if (index < (5+SoArealData::count+SoAintData::count)) {

            const int idx = index-(5+SoArealData::count);
            return particle.idata(idx);

          } else if (index < (5+SoArealData::count+SoAintData::count+m_rrData_count)) {

            const int idx = index-(5+SoArealData::count+SoAintData::count);
            return ptile_data.m_runtime_rdata[idx][i];
          }

          amrex::Abort("Invalid index for accessing particles data");
          return 0.;
        }

      private:
        const int m_rrData_count;
    };

    // The function to count the number of particles
    virtual
    void count () override;

  protected:
    // Set variables, components and indexes
    virtual
    void setup_variables () override;

    // The indexes
    amrex::Vector<int> m_components;

    // Setup the monitor
    virtual
    void initialize () override
    {
      Monitor::initialize();

      this->check_realbox_is_ok();
      m_direction = 0;
    }

    // Pointer to the particle container
    MFIXParticleContainer* m_pc;

    // Reference to the MFIXSolidsPhase class
    MFIXSolidsPhase* m_solids;

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override = 0;

    // Check the realbox has positive volume
    virtual
    void check_realbox_is_ok () const;
};

// General property class
class GeneralProperty final: public BaseMonitor {

  public:
    // Constructor
    GeneralProperty (MFIXParticleContainer* pc,
                     amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                     amrex::RealBox const& region)
      : BaseMonitor(pc, ebfactory, region)
    {}

    // Constructor
    GeneralProperty (const std::array<std::string,2> specs,
                     const std::string& pparse,
                     MFIXParticleContainer* pc,
                     amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                     MFIXSolidsPhase& solids,
                     MFIXRegions const& regions)
      : BaseMonitor(specs, pparse, pc, ebfactory, solids, regions)
    {}

    // The sum of particles property phi_p in the selected region
    amrex::Vector<amrex::Real> sum (const int lev,
                                    const amrex::Vector<int>& indexes);

    // The min of particles property phi_p in the selected region
    amrex::Vector<amrex::Real> min (const int lev,
                                    const amrex::Vector<int>& indexes);

    // The max of particles property phi_p in the selected region
    amrex::Vector<amrex::Real> max (const int lev,
                                    const amrex::Vector<int>& indexes);

    // The sum of particles property phi_p in the selected region
    amrex::Vector<amrex::Real> massweightedsum (const int lev,
                                                const amrex::Vector<int>& indexes);

  protected:
    // Setup the monitor
    virtual
    void initialize () override
    {
      BaseMonitor::initialize();
    }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override;
};

class AveragedProperty final: public BaseMonitor {

  public:
    // Constructor
    AveragedProperty (MFIXParticleContainer* pc,
                      amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                      amrex::RealBox const& region)
      : BaseMonitor(pc, ebfactory, region)
    {}

    // Constructor
    AveragedProperty (const std::array<std::string,2> specs,
                      const std::string& pparse,
                      MFIXParticleContainer* pc,
                      amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                      MFIXSolidsPhase& solids,
                      MFIXRegions const& regions)
      : BaseMonitor(specs, pparse, pc, ebfactory, solids, regions)
    {}

    // The average value of particles property phi_p in the selected region. For
    // DEM simulations the statistical weight of a particle is 1, such that the
    // sum of the weights is the total number of observations in the selected
    // region
    amrex::Vector<amrex::Real> average (const int lev,
                                        const amrex::Vector<int>& indexes);

    // The standard deviation of particles property phi_p in the selected region
    amrex::Vector<amrex::Real> stddev (const int lev,
                                       const amrex::Vector<int>& indexes);

    // Mass-weighted average value of particles property phi_p in the selected
    // region
    amrex::Vector<amrex::Real> mass_weighted_average (const int lev,
                                                      const amrex::Vector<int>& indexes);

    // Volume-weighted average value of particles property phi_p in the selected
    // region
    amrex::Vector<amrex::Real> volume_weighted_average (const int lev,
                                                        const amrex::Vector<int>& indexes);

  protected:
    // Setup the monitor
    virtual
    void initialize () override
    {
      BaseMonitor::initialize();
    }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override;
};

class FlowRate final: public BaseMonitor {

  public:
    // Constructor
    FlowRate (MFIXParticleContainer* pc,
              amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
              amrex::RealBox const& region,
              amrex::RealBox const* plane)
      : BaseMonitor(pc, ebfactory, region, plane)
    {
      this->check_plane_is_ok();
      this->set_direction();
      this->set_coordinate();
    }

    // Constructor
    FlowRate (const std::array<std::string,2> specs,
              const std::string& pparse,
              MFIXParticleContainer* pc,
              amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
              MFIXSolidsPhase& solids,
              MFIXRegions const& regions)
      : BaseMonitor(specs, pparse, pc, ebfactory, solids, regions)
    {}

    // The net flow rate of a general particle property phi_p is computed by
    // summing the properties of the set of particles projected to have crossed
    // the flow plane Gamma
    amrex::Vector<amrex::Real> flow_rate (const int lev,
                                          const amrex::Vector<int>& indexes,
                                          const amrex::Real dt);

    // The net mass-weighted flow rate is the sum of the general particles
    // property phi_p multiplied by the particles mass m_p of the set of
    // particles projected to have crossed the flow plane Gamma
    amrex::Vector<amrex::Real> mass_weighted_flow_rate (const int lev,
                                                        const amrex::Vector<int>& indexes,
                                                        const amrex::Real dt);

    // The net volume-weighted flow rate is the sum of the general particles
    // property phi_p multiplied by the particle volume V_p of the set of
    // particles projected to have crossed the flow plane Gamma
    amrex::Vector<amrex::Real> volume_weighted_flow_rate (const int lev,
                                                          const amrex::Vector<int>& indexes,
                                                          const amrex::Real dt);

    struct CrossesFlowPlane {

      AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
      int operator() (const amrex::Real& position,
                      const amrex::Real& velocity,
                      const amrex::Real& plane_coordinate,
                      const amrex::Real& dt) const;
    };

  protected:
    // Function for testing if particle crosses the flow plane
    // Setup the monitor
    virtual
    void initialize () override
    {
      BaseMonitor::initialize();

      m_plane = m_regions->getRegion(m_plane_name);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_plane != nullptr, "Region does not exist");

      this->check_plane_is_ok();
      this->set_direction();
      this->set_coordinate();
    }

    // The function to actually execute monitoring
    virtual
    void monitor (const amrex::Real& dt) override;

    // Check the box is good
    void check_plane_is_ok () const;

    // Function to set the area normal vector direction and check the box is 2D
    void set_direction ();

    // Function to set the plane coordinate along m_direction
    void set_coordinate ();

    // The plane coordinate along m_direction
    amrex::Real m_coordinate;
};

}

#endif
