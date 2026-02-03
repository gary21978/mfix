#ifndef MFIX_CHEMISTRY_MANAGEMENT_INTERFACE_H
#define MFIX_CHEMISTRY_MANAGEMENT_INTERFACE_H

#include <AMReX_REAL.H>
#include <AMReX_FabArray.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <integrator.H>

#include <mfix_pc.H>
#include <mfix_reactions.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_timer.H>


/**
 * \brief This class manages the whole chemistry in MFIX-Exa
 *
 * The ChemistryManagementInterface class contains all the routines needed in
 * MFIX-Exa to perform the chemistry update for the solid particles and compute
 * the chemistry transfer quantities for the fluid phase that will be used
 * during the fluid updating steps.
 */
class ChemistryManagementInterface
{
  public:
    /**
     * \brief Constructs an empty ChemistryManagementInterface
     *
     * Data can be defined at a later time using the define member function.
     */
    ChemistryManagementInterface ()
#ifdef OPTIMIZED_USR_RATES
      : rates_file_already_checked(0)
#endif
    {}

    /**
     * \brief Defines the ChemistryManagementInterface class
     *
     * \param reactions is the pointer to the reactions class
     * \param solids is the pointer to the solids class
     * \param fluid is the pointer to the fluid class
     * \param dem is the pointer to the DEM class
     * \param pic is the pointer to the PIC class
     * \param pc is the pointer to the NeighborParticleContainer class
     * \param timer is the pointer to the timer class
     */
    void update_pointers (MFIXReactions* reactions,
                          MFIXSolidsPhase* solids,
                          MFIXFluidPhase* fluid,
                          MFIXDEM* dem,
                          MFIXPIC* pic,
                          MFIXParticleContainer* pc,
                          MFIXTimer* timer,
                          const int report_mass_balance)
    {
      m_reactions = reactions;
      m_solids = solids;
      m_fluid = fluid;
      m_dem = dem;
      m_pic = pic;
      m_pc = pc;
      m_timer = timer;
      m_report_mass_balance = report_mass_balance;
    }

    /**
     * \brief This member function performs the whole chemistry update
     *
     * This function first computes all the information regarding which
     * particles are contained in which fluid cells, and then computes the
     * chemistry update
     *
     * \param geom is the reference to the problem's geometry
     * \param grids is the reference to the problem's grids
     * \param dmap is the reference to the problem's distribution map
     * \param particle_ebfactory is the reference to the problem's particles
     * EBFactory
     * \param txfr_out is the reference to the MultiFab where the fluid
     * chemistry transfer quantities will be stored
     * \param interp_ptr is the reference to the MultiFabs that contain all the
     * fluid information needed during the chemistry computation
     * \param thermodynamic_p_g_in is the fluid thermodynamic pressure
     * \param flags is the reference to the problem's EBCellFlagFab
     * \param volfrac is the reference to the fluid volume fraction MultiFab
     * \param time is the current simulation time
     * \param dt is the current simulation timestep that will be used for the
     * chemistry ODE integration
     * \param verbose is the flag for chemistry update verbosity
     */
    void calc_chemistry (amrex::Vector< amrex::Geometry > const& geom,
                         amrex::Vector< amrex::MultiFab* > const& chem_txfr_out,
                         amrex::Vector< amrex::MultiFab* > const& interp_ptr,
                         amrex::Real const& therm_p_in,
                         const amrex::Real time,
                         const amrex::Real dt,
                         const int verbose);

    /**
     * \brief Compute the information regarding which particles are contained in
     * which fluid cells.
     *
     * This member function computes the particles information for each
     * fluid cell by storing three data structures that contain the number of
     * particles in each fluid cell, a vector ordered by fluid cell that
     * contains (grouped by fluid cells) the particles IDs that are needed to
     * access the NeighborParticleContainer data, and another data struct which
     * provides the mapping between the fluid cells and the vector containing
     * the particles IDs. Indeed, as each fluid cell can have an unknown number
     * of particles, we need to map each fluid cell with its starting address in
     * the vector of particles IDs.
     *
     * \param time is the current simulation time
     * \param geom is the reference to the problem's geometry
     * \param interp_ptr is the reference to the MultiFabs that contain all the
     * fluid
     * \param parts_nb_in_cells is the number of particles in each fluid cell
     * \param parts_idxs_in_cells is the map between the fluid cells and the
     * corresponding first address in the vector containing the particles IDs
     * \param parts_idxs_vectorized is the monodimensional vector of particles
     * IDs ordered and grouped by fluid cell
     */
    void calc_particles_within_cells (const amrex::Real time,
                                      amrex::Vector<amrex::Geometry> const& geom,
                                      amrex::Vector< amrex::MultiFab* > const& interp_ptr,
                                      amrex::Vector< amrex::iMultiFab* > const& parts_nb_in_cells,
                                      amrex::Vector< amrex::iMultiFab* > const& parts_idxs_in_cells,
                                      amrex::Vector<std::map<MFIXParticleContainer::PairIndex,
                                                             amrex::Gpu::DeviceVector<int>>>& parts_idxs_vectorized);

    /**
     * \brief Compute the chemistry update
     *
     * This member function calls the actual chemistry update by passing all its
     * arguments plus the ODE integrator chosen by the user in the inputs file
     *
     * \param time is the current simulation time
     * \param dt is the current simulation timestep that will be used for the
     * chemistry ODE integration
     * \param geom is the reference to the problem's geometry
     * \param txfr_out is the reference to the MultiFab where the fluid
     * chemistry transfer quantities will be stored
     * \param interp_ptr is the reference to the MultiFabs that contain all the
     * fluid
     * \param therm_p_in is the fluid thermodynamic pressure
     * \param flags is the reference to the problem's EBCellFlagFab
     * \param volfrac is the reference to the fluid volume fraction MultiFab
     * \param parts_nb_in_cells is the number of particles in each fluid cell
     * \param parts_idxs_in_cells is the map between the fluid cells and the
     * corresponding first address in the vector containing the particles IDs
     * \param parts_idxs_vectorized is the monodimensional vector of particles
     * IDs ordered and grouped by fluid cell
     */
    void calc_chem_txfr (const amrex::Real time,
                         const amrex::Real dt,
                         amrex::Vector< amrex::Geometry > const& geom,
                         amrex::Vector< amrex::MultiFab* > const& txfr_out,
                         amrex::Vector< amrex::MultiFab* > const& interp_ptr,
                         amrex::Real const& therm_p_in,
                         amrex::Vector< amrex::iMultiFab* > const& parts_nb_in_cells,
                         amrex::Vector< amrex::iMultiFab* > const& parts_idxs_in_cells,
                         amrex::Vector< std::map< MFIXParticleContainer::PairIndex,
                                                  amrex::Gpu::DeviceVector<int> > >& parts_idxs_vectorized);

    /**
     * \brief Compute the chemistry update
     *
     * \param time is the current simulation time
     * \param dt is the current simulation timestep that will be used for the
     * chemistry ODE integration
     * \param geom is the reference to the problem's geometry
     * \param txfr_out is the reference to the MultiFab where the fluid
     * chemistry transfer quantities will be stored
     * \param interp_ptr is the reference to the MultiFabs that contain all the
     * fluid
     * \param therm_p_in is the fluid thermodynamic pressure
     * \param flags is the reference to the problem's EBCellFlagFab
     * \param volfrac is the reference to the fluid volume fraction MultiFab
     * \param parts_nb_in_cells is the number of particles in each fluid cell
     * \param parts_idxs_in_cells is the map between the fluid cells and the
     * corresponding first address in the vector containing the particles IDs
     * \param parts_idxs_vectorized is the monodimensional vector of particles
     * IDs ordered and grouped by fluid cell
     * \param integrator is the ODE integrator
     */
    template <class F1>
    void calc_chem_txfr (const amrex::Real time,
                         const amrex::Real dt,
                         amrex::Vector< amrex::Geometry > const& geom,
                         amrex::Vector< amrex::MultiFab* > const& txfr_out,
                         amrex::Vector< amrex::MultiFab* > const& interp_ptr,
                         amrex::Real const& therm_p_in,
                         amrex::Vector< amrex::iMultiFab* > const& parts_nb_in_cells,
                         amrex::Vector< amrex::iMultiFab* > const& parts_idxs_in_cells,
                         amrex::Vector< std::map< MFIXParticleContainer::PairIndex,
                                                  amrex::Gpu::DeviceVector<int> > >& parts_idxs_vectorized,
                         const F1& integrator);

    amrex::Gpu::DeviceVector<amrex::Real> const& fluid_mass_change () const
    { return m_fluid_mass_change; }

    amrex::Gpu::DeviceVector<amrex::Real> const& solids_mass_change () const
    { return m_solids_mass_change; }

  private:
#ifdef OPTIMIZED_USR_RATES
    int rates_file_already_checked;
#endif

    IntegratorInputs m_integrator_inputs;

    MFIXReactions* m_reactions;

    MFIXSolidsPhase* m_solids;

    MFIXFluidPhase* m_fluid;

    MFIXDEM* m_dem;

    MFIXPIC* m_pic;

    MFIXParticleContainer* m_pc;

    MFIXTimer* m_timer;

    int m_report_mass_balance;

    amrex::Gpu::DeviceVector<amrex::Real> m_fluid_mass_change;
    amrex::Gpu::DeviceVector<amrex::Real> m_solids_mass_change;
};


#endif
