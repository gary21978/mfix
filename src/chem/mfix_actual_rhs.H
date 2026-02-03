#ifndef MFIX_ACTUAL_RHS_H
#define MFIX_ACTUAL_RHS_H

#include <type_traits>

#include <AMReX_REAL.H>

#include <mfix_chem_data_structs.H>
#include <mfix_usr_reactions_rates_K.H>


/**
 * \brief This class computes the right-hand side of the chemistry ODE to be
 * integrated.
 *
 * In each fluid cell and for each particle in the fluid cell, the ODE system we
 * integrate is:
 * \f[
 * \left\{ \begin{array}{l}
 *   \displaystyle\frac{D \left(\rho_g X_{gk} \right)}{D t} =
 *       \displaystyle\frac{\mathcal{G}^{\rho}_{p \, \rightarrow gk}}{\varepsilon_g},
 *       \quad \forall k = 1,\ldots,N_g \\
 *   \displaystyle\frac{D \left(\rho_g h_g \right)}{D t} =
 *       \displaystyle\frac{\mathcal{G}^{h}_{p \, \rightarrow g}}{\varepsilon_g} \\
 *   \displaystyle\frac{d \left(m_p X_{pk}\right)}{dt} =
 *       G^{m}_{g \, \rightarrow pk}, \quad \forall k = 1,\ldots,N_s \\
 *   \displaystyle\frac{d\left(m_p h_p\right)}{dt} =
 *       G^{H}_{g \, \rightarrow p}
 * \end{array} \right.
 * \f]
 * where \f$ \varepsilon_g \f$, \f$ \rho_g \f$, \f$ h_g \f$, and \f$ X_{gk} \f$
 * are respectively the fluid volume fraction, density, specific enthalpy, and
 * mass fractions; \f$ m_p \f$, \f$ h_p \f$, and \f$ X_{pk} \f$ are respectively
 * the particle mass, specific enthalpy, and mass fractions. \f$ N_g \f$ and \f$
 * N_s \f$ are respectively the number of species for the fluid and the solids
 * phases. \f$ \mathcal{G}^{\rho}_{p \, \rightarrow gk} \f$ is the interphase
 * mass transfer rate from the particle to the k-th fluid species; \f$
 * \mathcal{G}^{h}_{p \, \rightarrow g} \f$ is the interphase specific enthalpy
 * transfer rate to the fluid phase due to the chemical reactions; \f$ G^{m}_{g
 * \, \rightarrow pk} \f$ is the mass transfer rate from the fluid phase to the
 * particle k-th species; \f$ G^{H}_{g \, \rightarrow p} \f$ is the interphase
 * enthalpy transfer rate to the particle due to the chemical reactions.
 *
 * This class gets the current state for the fluid phase and the particle,
 * computes the ODE right-hand side, and store it into the corresponding memory
 * areas.
 */
class MFIXActualRHS
{
  public:
    using ChemFluidIndexes = ChemMemoryIndexes::ChemFluidIndexes;
    using ChemSolidsIndexes = ChemMemoryIndexes::ChemSolidsIndexes;

    /**
     * \brief Constructs an empty MFIXActualRHS
     *
     * Data can be defined at a later time using the define member function.
     */
    AMREX_GPU_HOST_DEVICE
    MFIXActualRHS ()
      : m_rate_type(ReactionType::Undefined)
      , m_fluid_data(nullptr)
      , m_solids_data(nullptr)
      , m_reactions_data(nullptr)
      , m_indexes(nullptr)
    {}

    /**
     * \brief Defines the MFIXActualRHS class
     *
     * This member function assigns the class members for a correct computation
     * of the ODE right-hand side
     *
     * \param rate_type is the flag to select which type of reactions to use for
     * RHS computation
     * \param fluid_data is the class containing fluid variables and information
     * needed for RHS computation
     * \param solids_data is the class containing current particle's variables
     * and information needed for RHS computation
     * \param reactions_data is the class containing reactions information
     * needed for RHS computation
     * \param indexes is the class containing indexing information for properly
     * access fluid and solids memory speces
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup (const int rate_type,
                const int solve_fluid,
                const int solve_solids,
                MFIXFluidData* fluid_data,
                MFIXSolidsData* solids_data,
                MFIXReactionsData* reactions_data,
                const ChemMemoryIndexes& indexes)
    {
      // Assign the class members
      m_rate_type = rate_type;
      m_solve_fluid = solve_fluid;
      m_solve_solids = solve_solids;
      m_fluid_data = fluid_data;
      m_solids_data = solids_data;
      m_reactions_data = reactions_data;
      m_indexes = &indexes;
    }

    /**
     * \brief Compute the ODE right-hand side
     *
     * Function call operator for the MFIXActualRHS class which computes the ODE
     * right-hand side
     *
     * \param y is the ODE solution
     * \param y_rhs is the ODE right-hand side
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (ChemVector<amrex::Real>& y,
                      ChemVector<amrex::Real>& y_rhs) const
    {
      // Reset the right hand side
      y_rhs.fill(0.);

      const amrex::Real& ep_g = m_solve_fluid ? m_fluid_data->volume_fraction() : 0.;
      const amrex::Real& ro_g = m_solve_fluid ? m_fluid_data->density() : 0.;

      // Get the number of reactions
      const int nreactions = m_reactions_data->params.m_nreactions;

      // if the reactions are purely eulerian
      if (isEulerian() && m_solve_fluid) {

        // Set fluid variables to be equal to argument variables. This is
        // important for iterative ODE integrators as they call the RHS also on
        // temporary data and not just on the problem unknowns. The fluid
        // variables set by set_y_ptr and set_y_rhs_ptr are used in the body of
        // this call operator function
        m_fluid_data->set_y_ptr(y.ptr(m_indexes->fluid_first_idx()));
        m_fluid_data->set_y_rhs_ptr(y_rhs.ptr(m_indexes->fluid_first_idx()));

        // Compute reactions rates
        m_eulerian_rates(*m_reactions_data, *m_fluid_data);

        amrex::Real G_rho_g(0);

        amrex::Real const Tg = m_fluid_data->temperature();

        amrex::Real cp_g(0.);

        // Loop over fluid species
        for (int n_g(0); n_g < m_fluid_data->params.nspecies(); ++n_g) {

          amrex::Real G_rho_gk(0);

          amrex::Real const MW_gk = m_fluid_data->molar_mass(n_g);

          amrex::Real const h_gk = m_fluid_data->m_props.enthalpy(n_g, Tg);

          cp_g += m_fluid_data->mass_fraction(n_g) *
                  m_fluid_data->m_props.specificHeat(n_g, Tg);

          // Loop over the reactions to compute each contribution to the ODE
          // right-hand side
          for (int q(0); q < nreactions; q++) {

            // Total rate of change for n_g-th fluid species mass (kg/sec/m3)
            // Multiplying by scale converts to "per volume"
            G_rho_gk += m_fluid_data->params.get_eulerian_stoich_coeff(n_g,q) *
                        m_reactions_data->m_rates[q] * MW_gk;
          }

          G_rho_g += G_rho_gk;

          // Update fluid n_g-th species density RHS with current reaction
          // contribution
          m_fluid_data->m_rhs[ChemFluidIndexes::mass_fractions()+n_g] = G_rho_gk;

          m_fluid_data->m_rhs[ChemFluidIndexes::temperature()] -= G_rho_gk * h_gk;
        }

        //
        m_fluid_data->m_rhs[ChemFluidIndexes::density()] = (1./ep_g) * G_rho_g;

        //
        for (int n_g(0); n_g < m_fluid_data->params.nspecies(); ++n_g) {

          amrex::Real G_rho_gk = m_fluid_data->m_rhs[ChemFluidIndexes::mass_fractions()+n_g];
          amrex::Real X_gk = m_fluid_data->mass_fraction(n_g);

          m_fluid_data->m_rhs[ChemFluidIndexes::mass_fractions()+n_g] =
            (1./(ep_g*ro_g)) * (G_rho_gk - X_gk * G_rho_g);
        }

        m_fluid_data->m_rhs[ChemFluidIndexes::temperature()] /= (ep_g*ro_g*cp_g);

        m_fluid_data->m_rhs[ChemFluidIndexes::specific_enthalpy()] = 0.;

        m_fluid_data->reset_y_ptr();
        m_fluid_data->reset_y_rhs_ptr();

      } else if (isLagrangian() && m_solve_solids) { // if rate type flag is lagrangian

        const amrex::Real mass_p = m_solids_data->mass();

        if (m_solve_fluid) {
          // Set fluid variables to be equal to argument variables. This is
          // important for iterative ODE integrators as they call the RHS also
          // on temporary data and not just on the problem unknowns. The fluid
          // variables set by set_y_ptr and set_y_rhs_ptr are used in the body of
          // this call operator function
          m_fluid_data->set_y_ptr(y.ptr(m_indexes->fluid_first_idx()));
          m_fluid_data->set_y_rhs_ptr(y_rhs.ptr(m_indexes->fluid_first_idx()));
        }

        // Set solids variables to be equal to argument variables. This is
        // important for iterative ODE integrators as they call the RHS also on
        // temporary data and not just on the problem unknowns. The solids
        // variables set by set_y_ptr and set_y_rhs_ptr are used in the body of
        // this call operator function
        m_solids_data->set_y_ptr(y.ptr(m_indexes->solids_first_idx()));
        m_solids_data->set_y_rhs_ptr(y_rhs.ptr(m_indexes->solids_first_idx()));

        // Compute reactions rates
        m_lagrangian_rates(*m_reactions_data, *m_fluid_data, *m_solids_data);

        // Local auxiliary variables for fluid and current particle's
        // temperature
        amrex::Real const Ts(m_solids_data->temperature());
        amrex::Real const Tg = m_solve_fluid?
          m_fluid_data->temperature() : std::numeric_limits<amrex::Real>::infinity();

        // Local variable to store heat of reaction
        amrex::Real HoR(0.);

        amrex::Real cp_g(0.);
        const amrex::Real h_g = m_fluid_data->specific_enthalpy();

        amrex::Real G_rho_g(0.);

        // Loop over fluid phase species
        for (int n_g(0); n_g < m_fluid_data->params.nspecies() && m_solve_fluid; ++n_g) {

          amrex::Real G_m_gk(0);

          amrex::Real const MW_gk = m_fluid_data->molar_mass(n_g);

          amrex::Real const h_gk_Ts = m_fluid_data->m_props.enthalpy(n_g, Ts);
          amrex::Real const h_gk_Tg = m_fluid_data->m_props.enthalpy(n_g, Tg);

          cp_g += m_fluid_data->mass_fraction(n_g) * m_fluid_data->m_props.specificHeat(n_g, Tg);

          // Loop over the reactions to compute each contribution to the ODE
          // right-hand side
          for (int q(0); q < nreactions; q++) {

            // Total rate of change for n_g-th fluid species mass (kg/sec)
            amrex::Real G_m_gk_q = m_fluid_data->params.get_lagrangian_stoich_coeff(n_g, q) *
                                   m_reactions_data->m_rates[q] * MW_gk;

            G_m_gk += G_m_gk_q;

            if (m_reactions_data->has_fluid_and_solids(q)) {
              // Update the heat of reaction for his reaction
              HoR += (G_m_gk_q > 0.) ? G_m_gk_q*h_gk_Ts : G_m_gk_q*h_gk_Tg;
            }
          }

          G_rho_g += toEulerian(G_m_gk);

          // Update the right-hand side of the n_g-th fluid species density by
          // adding this reaction's contribution converted to "per volume"
          m_fluid_data->m_rhs[ChemFluidIndexes::mass_fractions()+n_g] = toEulerian(G_m_gk);

          //
          m_fluid_data->m_rhs[ChemFluidIndexes::temperature()] -= toEulerian(G_m_gk) * h_gk_Tg;

          m_fluid_data->m_rhs[ChemFluidIndexes::specific_enthalpy()] -= h_g * toEulerian(G_m_gk);
        }

        //
        m_fluid_data->m_rhs[ChemFluidIndexes::density()] = (1./ep_g) * G_rho_g;

        //
        for (int n_g(0); n_g < m_fluid_data->params.nspecies(); ++n_g) {

          amrex::Real G_rho_gk = m_fluid_data->m_rhs[ChemFluidIndexes::mass_fractions()+n_g];
          amrex::Real X_gk = m_fluid_data->mass_fraction(n_g);

          m_fluid_data->m_rhs[ChemFluidIndexes::mass_fractions()+n_g] =
            (1./(ep_g*ro_g)) * (G_rho_gk - X_gk * G_rho_g);
        }

        amrex::Real cp_p(0.);
        const amrex::Real h_p = m_solids_data->specific_enthalpy();

        amrex::Real G_m_s(0);

        // Loop over solids phase species
        for (int n_s(0); n_s < m_solids_data->params.nspecies(); ++n_s) {

          amrex::Real G_m_sn(0);

          amrex::Real const MW_sn = m_solids_data->molar_mass(n_s);

          amrex::Real const h_sn_Tg = m_solids_data->m_props.enthalpy(n_s, Tg);
          amrex::Real const h_sn_Ts = m_solids_data->m_props.enthalpy(n_s, Ts);

          cp_p += m_solids_data->mass_fraction(n_s) * m_solids_data->m_props.specificHeat(n_s, Ts);

          for (int q(0); q < nreactions; q++) {

            amrex::Real const& rate = m_reactions_data->m_rates[q];

            // Total rate of change for the n_s-th solids species for current
            // particle (kg/sec)
            amrex::Real const G_m_sn_q = rate * MW_sn
              * m_solids_data->params.get_lagrangian_stoich_coeff(n_s,q);

            G_m_sn += G_m_sn_q;

            if (m_solve_fluid && (m_reactions_data->has_fluid_and_solids(q))) {
              // Update the heat of reaction for this reaction
              HoR += G_m_sn_q * (h_sn_Tg - h_sn_Ts);
            }
          }

          G_m_s += G_m_sn;

          // Update the right-hand side of the n_s-th solids species mass for the
          // current particle by adding this reaction's contribution
          m_solids_data->m_rhs[ChemSolidsIndexes::mass_fractions()+n_s] = G_m_sn;

          //
          m_solids_data->m_rhs[ChemSolidsIndexes::temperature()] -= G_m_sn * h_sn_Ts;

          m_solids_data->m_rhs[ChemSolidsIndexes::specific_enthalpy()] -= h_p * G_m_sn;
        }

        // if constant particle volume
        if (m_solids_data->constant_volume()) {
          m_solids_data->m_rhs[ChemSolidsIndexes::radius()] = 0;
          m_solids_data->m_rhs[ChemSolidsIndexes::density()] = G_m_s / m_solids_data->volume();
        } else { // constant density
          amrex::Real radius_square = m_solids_data->radius()*m_solids_data->radius();
          amrex::Real coeff = 4*M_PI*radius_square*m_solids_data->density();
          m_solids_data->m_rhs[ChemSolidsIndexes::radius()] = G_m_s / coeff;
          m_solids_data->m_rhs[ChemSolidsIndexes::density()] = 0;
        }

        //
        for (int n_s(0); n_s < m_solids_data->params.nspecies(); ++n_s) {

          amrex::Real G_m_sn = m_solids_data->m_rhs[ChemSolidsIndexes::mass_fractions()+n_s];
          amrex::Real X_sn = m_solids_data->mass_fraction(n_s);

          m_solids_data->m_rhs[ChemSolidsIndexes::mass_fractions()+n_s] =
            (1./mass_p) * (G_m_sn - X_sn * G_m_s);
        }

        // Update the right-hand side of the fluid and solids enthalpy equations
        // in the ODE to be integrated. Note that in the fluid equation the heat
        // of reaction is converted to "per volume"
        if (m_solve_fluid) {
          m_fluid_data->m_rhs[ChemFluidIndexes::specific_enthalpy()] += toEulerian(HoR);
          m_fluid_data->m_rhs[ChemFluidIndexes::specific_enthalpy()] /= (ep_g*ro_g);

          m_fluid_data->m_rhs[ChemFluidIndexes::temperature()] += toEulerian(HoR);
          m_fluid_data->m_rhs[ChemFluidIndexes::temperature()] /= (ep_g*ro_g*cp_g);
        }

        m_solids_data->m_rhs[ChemSolidsIndexes::specific_enthalpy()] -= HoR;
        m_solids_data->m_rhs[ChemSolidsIndexes::specific_enthalpy()] /= (mass_p);

        m_solids_data->m_rhs[ChemSolidsIndexes::temperature()] -= HoR;
        m_solids_data->m_rhs[ChemSolidsIndexes::temperature()] /= (mass_p*cp_p);

        if (m_solve_fluid) {
          // Reset the fluid variables to be equal to te problem unknowns in
          // case we have changed them at the beginning of this function to be
          // something different
          m_fluid_data->reset_y_ptr();
          m_fluid_data->reset_y_rhs_ptr();
        }

        // Reset the solids variables to be equal to te problem unknowns in case
        // we have changed them at the beginning of this function to be
        // something different
        m_solids_data->reset_y_ptr();
        m_solids_data->reset_y_rhs_ptr();

      } // end if statement on eulerian vs eagrangian rates
    } // end loop q-th reaction

  private:
    int m_rate_type;

    int m_solve_fluid;
    int m_solve_solids;

    MFIXFluidData* m_fluid_data;
    MFIXSolidsData* m_solids_data;
    MFIXReactionsData* m_reactions_data;
    const ChemMemoryIndexes* m_indexes;

    LagrangianReactionRates m_lagrangian_rates;
    EulerianReactionRates m_eulerian_rates;

    /**
     * Returns true if the current rates type is eulerian, false otherwise
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int isEulerian () const noexcept
    { return m_rate_type == ReactionType::Eulerian ? 1 : 0; }

    /**
     * Returns true if the current rates type is lagrangian, false otherwise
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int isLagrangian () const noexcept
    { return m_rate_type == ReactionType::Eulerian ? 0 : 1; }

    /**
     * \brief Convert lagrangian quantities to eulerian
     *
     * Member function that transforms lagrangian quatities to eulerian
     * quantities by multiplying by the statistical weight (in case of PIC), and
     * by dividing by the effective cell volume
     *
     * \param a_rate is the quantity to be converted
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real toEulerian ( amrex::Real const a_rate  ) const noexcept
    {
      amrex::Real scaled_rate(a_rate);
      scaled_rate *= m_solids_data->statistical_weight();
      scaled_rate /= (m_reactions_data->cell_volume()
                      * m_reactions_data->vfrac());
      return scaled_rate;
    }
};


#endif
