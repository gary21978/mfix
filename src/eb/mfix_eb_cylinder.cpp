#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <algorithm>
#include <mfix_eb.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple cylinder EB.                                     *
 *                                                                              *
 * Comments: The EB for the cylinder is open at the top and bottom; if the      *
 * boundary condition is mass inflow then we modify the level set so that       *
 * the particles see it as a solid boundary.                                    *
 *                                                                              *
 *******************************************************************************/
void MFIXEB::make_eb_cylinder (Vector<Geometry> const& a_geom)
{
    ParmParse pp("cylinder");

    /****************************************************************************
     * Get cylinder information from inputs file.                               *
     ***************************************************************************/
    bool inside       = true;

    Real radius    = 0.0002;
    Real height    = -1.;

    int direction  = 0;
    Vector<Real> centervec(3);
    Real rotation = 0.;
    int rotation_axis = 0;

    pp.query("internal_flow", inside);

    pp.query("radius",       radius);
    pp.query("height",       height);
    pp.query("direction",    direction);
    pp.getarr("center",      centervec,  0, 3);
    pp.query("rotation",     rotation);
    pp.query("rotation_axis", rotation_axis);
    Array<Real,3> center={centervec[0], centervec[1], centervec[2]};


    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    amrex::Print() << " " << std::endl;
    amrex::Print() << " Internal Flow: " << inside << std::endl;
    amrex::Print() << " Radius:    " << radius    << std::endl;
    amrex::Print() << " Height:    " << height    << std::endl;
    amrex::Print() << " Direction: " << direction << std::endl;
    amrex::Print() << " Rotation angle: " << rotation << std::endl;
    amrex::Print() << " Rotation axis: " << rotation_axis << std::endl;
    amrex::Print() << " Center:    " << center[0] << ", "
                   << center[1] << ", "
                   << center[2] << std::endl;

    if(rotation == 0.0) {
       // Build the Cylinder geometry first representing the curved walls (this is
       // always present regardless of user input).
       EB2::CylinderIF my_cyl(radius, height, direction, center, inside);

       auto gshop_cyl = EB2::makeShop(my_cyl);

       build_levels(a_geom, gshop_cyl);
    }
    else {
       AMREX_ALWAYS_ASSERT_WITH_MESSAGE(height == -1, "You cannot rotate cylinders of finite length.");

       rotation = (rotation/180.) * M_PI;

       Array<Real,3> origin = {0.0, 0.0, 0.0};
       EB2::CylinderIF my_cyl(radius, height, direction, origin, inside);
       auto my_cyl_rotated = EB2::rotate(my_cyl, rotation, rotation_axis);
       auto my_cyl_final = EB2::translate(my_cyl_rotated,  center);

       auto gshop_cyl = EB2::makeShop(my_cyl_final);

       build_levels(a_geom, gshop_cyl);
    }

    Print() << "Done making the fluid eb levels ..." << std::endl;
}
