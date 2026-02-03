#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Rotation.H>

#include <algorithm>
#include <AMReX_EB_utils.H>
#include <mfix_eb.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple hopper EB.                                       *
 *                                                                              *
 * Comments: The hopper is constructed by "lathing" a union of two planes       *
 * (representing the hopper wall and outlet funnel) around a central axis.      *
 *                                                                              *
 *******************************************************************************/

void
MFIXEB::make_eb_hopper (Vector<Geometry> const& a_geom)
{
    ParmParse pp("hopper");

    /****************************************************************************
     * Get hopper information from inputs file.                                 *
     ***************************************************************************/

    Real orifice_radius  = 0.00025;
    Real funnel_radius   = 0.00050;
    Real funnel_height   = 0.00105;

    int direction = 0;
    Vector<Real> centervec(3);

    pp.query("orifice_radius", orifice_radius);
    pp.query("funnel_radius",  funnel_radius);
    pp.query("funnel_height",  funnel_height);

    pp.query("direction", direction);
    pp.getarr("center",   centervec,  0, 3);
    Array<Real,3> center={centervec[0],centervec[1],centervec[2]};

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    amrex::Print() << " " << std::endl;
    amrex::Print() << " Orifice Radius: " << orifice_radius << std::endl;
    amrex::Print() << " Funnel Radius:  " << funnel_radius  << std::endl;
    amrex::Print() << " Funnel Height:  " << funnel_height  << std::endl;
    amrex::Print() << " Direction: "      << direction      << std::endl;
    amrex::Print() << " Center:    "      << center[0]      << ", "
                   << center[1]      << ", "
                   << center[2]      << std::endl;


    // Construct the plane representing the funnel edge
    Array<Real,3> point{0.0, 0.0, 0.0};
    point[0] = orifice_radius;

    Array<Real,3> normal{0.0, 0.0, 0.0};
    normal[0] =  funnel_height;
    normal[1] = orifice_radius - funnel_radius;

    amrex::Print() << " " << std::endl;
    amrex::Print() << " Funnel plane: " << std::endl;
    amrex::Print() << "   Point:  " << point[0]  << ", "
                   << point[1]  << ", "
                   << point[2]  << std::endl;

    amrex::Print() << "   Normal: " << normal[0] << ", "
                   << normal[1] << ", "
                   << normal[2] << std::endl;

    EB2::PlaneIF funnel(point, normal);


    // Construct the plane representing the hopper-body's edge
    point = {funnel_radius, 0., 0.};
    normal = {1.0, 0.0, 0.0};

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Hopper body plane: " << std::endl;
    amrex::Print() << "   Point:  " << point[0]  << ", "
                   << point[1]  << ", "
                   << point[2]  << std::endl;

    amrex::Print() << "   Normal: " << normal[0] << ", "
                   << normal[1] << ", "
                   << normal[2] << std::endl;

    EB2::PlaneIF hbody(point, normal);

    // Basic hopper shape (by lathing the union of the body and hopper edges)
    auto hopper0 = EB2::lathe(EB2::makeUnion(hbody, funnel));

    // Align to correct axis
    auto hopper1 = EB2::rotate(hopper0, 2.*std::atan(1.0), 1);

    // Translate to correct location
    auto my_hopper = EB2::translate(hopper1, center);

    // Construct EB2 Index Space
    auto gshop = EB2::makeShop(my_hopper);

    build_levels(a_geom, gshop);

    Print() << "Done making the fluid eb levels ..." << std::endl;
}
