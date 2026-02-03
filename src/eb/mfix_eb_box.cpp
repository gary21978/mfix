#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>

#include <algorithm>
#include <AMReX_EB_utils.H>
#include <mfix_eb.H>
#include <mfix_eb_if.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple rectangular box with EB walls.                   *
 *                                                                              *
 *******************************************************************************/
void MFIXEB::make_eb_box (Vector<Geometry> const& a_geom)
{
    ParmParse pp_box("box");

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    // set up ebfactory

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Now making the ebfactories ..." << std::endl;

    if ( a_geom[0].isAllPeriodic() )
    {
        make_eb_regular(a_geom);
    }
    else
    {
        /************************************************************************
         *                                                                      *
         * Define Box geometry:                                                 *
         *        -> box.{Lo,Hi} vector storing box lo/hi                       *
         *        -> box.offset  vector storing box offset                      *
         * NOTE: walls are placed _outside_ domain for periodic directions.     *
         *                                                                      *
         ***********************************************************************/

        bool inside       = true;

        Vector<Real> boxLo(3), boxHi(3);
        Real offset    = 1.0e-15;

        for (int i = 0; i < 3; i++) {

            boxLo[i] = a_geom[0].ProbLo(i);
            boxHi[i] = a_geom[0].ProbHi(i);
        }

        pp_box.queryarr("Lo", boxLo,  0, 3);
        pp_box.queryarr("Hi", boxHi,  0, 3);

        pp_box.query("internal_flow", inside);
        pp_box.query("offset", offset);

        Real xlo = boxLo[0] + (inside ? offset : 0.);
        Real xhi = boxHi[0] - (inside ? offset : 0.);

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if (a_geom[0].isPeriodic(0) ) {

          if (boxLo[0] - offset <= a_geom[0].ProbLo(0))
          { xlo = 2.0*a_geom[0].ProbLo(0) - a_geom[0].ProbHi(0); }

          if (boxHi[0] + offset >= a_geom[0].ProbHi(0))
          { xhi = 2.0*a_geom[0].ProbHi(0) - a_geom[0].ProbLo(0); }
        }

        Real ylo = boxLo[1] + (inside ? offset : 0.);
        Real yhi = boxHi[1] - (inside ? offset : 0.);

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if (a_geom[0].isPeriodic(1)) {

          if (boxLo[1] - offset <= a_geom[0].ProbLo(1))
          { ylo = 2.0*a_geom[0].ProbLo(1) - a_geom[0].ProbHi(1); }

          if (boxHi[1] + offset >= a_geom[0].ProbHi(1))
          {  yhi = 2.0*a_geom[0].ProbHi(1) - a_geom[0].ProbLo(1); }
        }

        Real zlo = boxLo[2] + (inside ? offset : 0.);
        Real zhi = boxHi[2] - (inside ? offset : 0.);

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if (a_geom[0].isPeriodic(2)) {

          if (boxLo[2] - offset <= a_geom[0].ProbLo(1))
          { zlo = 2.0*a_geom[0].ProbLo(2) - a_geom[0].ProbHi(2); }

          if (boxHi[2] + offset >= a_geom[0].ProbHi(2))
          { zhi = 2.0*a_geom[0].ProbHi(2) - a_geom[0].ProbLo(2); }
        }

        Array<Real,3>  point_lox{ xlo, 0.0, 0.0};
        Array<Real,3> normal_lox{-1.0, 0.0, 0.0};
        Array<Real,3>  point_hix{ xhi, 0.0, 0.0};
        Array<Real,3> normal_hix{ 1.0, 0.0, 0.0};

        Array<Real,3>  point_loy{0.0, ylo, 0.0};
        Array<Real,3> normal_loy{0.0,-1.0, 0.0};
        Array<Real,3>  point_hiy{0.0, yhi, 0.0};
        Array<Real,3> normal_hiy{0.0, 1.0, 0.0};

        Array<Real,3>  point_loz{0.0, 0.0, zlo};
        Array<Real,3> normal_loz{0.0, 0.0,-1.0};
        Array<Real,3>  point_hiz{0.0, 0.0, zhi};
        Array<Real,3> normal_hiz{0.0, 0.0, 1.0};

        EB2::PlaneIF plane_lox(point_lox,normal_lox);
        EB2::PlaneIF plane_hix(point_hix,normal_hix);

        EB2::PlaneIF plane_loy(point_loy,normal_loy);
        EB2::PlaneIF plane_hiy(point_hiy,normal_hiy);

        EB2::PlaneIF plane_loz(point_loz,normal_loz);
        EB2::PlaneIF plane_hiz(point_hiz,normal_hiz);

        auto my_box = EB2::makeUnion(plane_lox, plane_hix, plane_loy,
                                     plane_hiy, plane_loz, plane_hiz );


        if( inside ) {

          auto gshop = EB2::makeShop(my_box);
          build_levels(a_geom, gshop);

        } else {

          auto gshop = EB2::makeShop(EB2::makeComplement(my_box));
          build_levels(a_geom, gshop);

        }

        Print() << "Done making the eb levels ..." << std::endl;
        Print() << " " << std::endl;
    }
}
