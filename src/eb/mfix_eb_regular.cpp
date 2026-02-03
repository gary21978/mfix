#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <mfix_eb.H>
#include <mfix_eb_if.H>

#include <AMReX_EB_utils.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Placeholder: create a simulation box _without_ EB walls.                     *
 *                                                                              *
 *******************************************************************************/
void
MFIXEB::make_eb_regular (Vector<Geometry> const& a_geom)
{
    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    // set up ebfactory

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Now making the ebfactories ..." << std::endl;

    // If filling level-set: this is used to store the implicit function (due to
    // any walls defined in mfix.dat). It is filled while after EB2::Build.
    // NOTE: this pointer is undefined if and of:
    //     * ! m_dem.solve()
    //     * levelset_restart
    // are true

    //MultiFab* mf_impfunc;

    // if (fluid.solve())
    EB2::AllRegularIF my_regular;
    auto gshop = EB2::makeShop(my_regular);

    build_levels(a_geom, gshop);
}
