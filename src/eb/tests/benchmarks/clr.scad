$fn = 50;
$fn = 50;

module air_reactor_and_riser(reactor_radius,
                             reactor_height,
                             riser_radius,
                             riser_height,
                             c2c_height) {
    overlap = 1e-4*reactor_height;
    cylinder(r=riser_radius, h=2*riser_height, center = true);
    cylinder(r=reactor_radius, h=2*(reactor_height + overlap), center = true);
    translate([0, 0, reactor_height]) {
        cylinder(r1=reactor_radius, r2=riser_radius, h=c2c_height);
    }
}

module crossover(crossover_radius,
                 crossover_height,
                 crossover_length) {

    translate([crossover_height, 0, 0]) {
        cylinder(r=crossover_radius, h=crossover_length);
    }
}

module lvalve(lvalve_radius,
              lvalve_height,
              lvalve_length) {

    translate([lvalve_height, 0, 0]) {
        cylinder(r=lvalve_radius, h=lvalve_length);
    }
}

module cyclone_wdipleg(cyclone_radius,
                       cyclone_height,
                       cyclone_dipleg_radius,
                       cyclone_dipleg_height,
                       cyclone_top,
                       cyclone_dipleg_bottom) {

    translate([0, 0, cyclone_dipleg_bottom]) {
        cylinder(h=2*cyclone_top, r=cyclone_dipleg_radius);
    }
    translate([0, 0, cyclone_top - cyclone_height - cyclone_dipleg_height]) {
        cylinder(h=cyclone_dipleg_height, r1=cyclone_dipleg_radius, r2=cyclone_radius);
        translate([0, 0, cyclone_dipleg_height]) {
            cylinder(h=cyclone_height, r=cyclone_radius);
        }
    }
}

module loopseal(loopseal_radius,
                loopseal_height,
                loopseal_bottom) {

    translate([0, 0, loopseal_bottom]) {
        cylinder(h=loopseal_height, r=loopseal_radius);
    }
}

module fuelreactor(fuelreactor_radius,
                   fuelreactor_height,
                   fuelreactor_bottom) {

    translate([0, 0, fuelreactor_bottom]) {
        cylinder(h=fuelreactor_height, r=fuelreactor_radius);
    }
}

module clr(reactor_radius,
           reactor_height,
           riser_radius,
           riser_height,
           c2c_height,
           offset,
           crossover_radius,
           crossover_height,
           crossover_length,
           cyclone_radius,
           cyclone_height,
           cyclone_dipleg_radius,
           cyclone_dipleg_height,
           cyclone_top,
           cyclone_dipleg_bottom,
           loopseal_radius,
           loopseal_height,
           loopseal_bottom,
           loopseal_offset,
           fuelreactor_radius,
           fuelreactor_height,
           fuelreactor_bottom,
           lvalve_radius,
           lvalve_height,
           lvalve_length) {

    translate([0, offset, offset]) {

        // Air reactor & riser
        rotate([0, 90, 0]) {
            air_reactor_and_riser(reactor_radius=reactor_radius,
                                  reactor_height=reactor_height,
                                  riser_radius=riser_radius,
                                  riser_height=riser_height,
                                  c2c_height=c2c_height);
        }

        // Crossover
        crossover(crossover_radius=crossover_radius,
                  crossover_height=crossover_height,
                  crossover_length=crossover_length);

        // Cyclone, loopseal and fuel reactor
        translate([0, 0, crossover_length]) {
            rotate([0, 90, 0]) {

                // Cylone with dipleg
                cyclone_wdipleg(cyclone_radius=cyclone_radius,
                                cyclone_height=cyclone_height,
                                cyclone_dipleg_radius=cyclone_dipleg_radius,
                                cyclone_dipleg_height=cyclone_dipleg_height,
                                cyclone_top=cyclone_top,
                                cyclone_dipleg_bottom=cyclone_dipleg_bottom);

                // Fuel reactor
                fuelreactor(fuelreactor_radius=fuelreactor_radius,
                            fuelreactor_height=fuelreactor_height,
                            fuelreactor_bottom=fuelreactor_bottom);
            }
        }

        // Loopseal
        translate([0, 0, crossover_length + loopseal_offset]) {
            rotate([0, 90, 0]) {
                // Loopseal
                loopseal(loopseal_radius=loopseal_radius,
                         loopseal_height=loopseal_height,
                         loopseal_bottom=loopseal_bottom);
            }
        }

        // Lvalve
        lvalve(lvalve_radius=lvalve_radius,
               lvalve_height=lvalve_height,
               lvalve_length=lvalve_length);
    }
}

module exa_clr() {
    clr(reactor_radius=0.1,
        reactor_height=0.75,
        riser_radius=0.05,
        riser_height=3.75,
        c2c_height=0.1,
        offset=0.16,
        crossover_radius=0.05,
        crossover_height=3.6,
        crossover_length=0.636,
        cyclone_radius=0.1,
        cyclone_height=0.35,
        cyclone_dipleg_radius=0.04,
        cyclone_dipleg_height=0.2,
        cyclone_top=3.8,
        cyclone_dipleg_bottom=0.15,
        loopseal_radius=0.12,
        loopseal_height=0.5,
        loopseal_bottom=2.24,
        loopseal_offset=0.004,
        fuelreactor_radius=0.12,
        fuelreactor_height=1.54,
        fuelreactor_bottom=0.52,
        lvalve_radius=0.05,
        lvalve_height=0.2,
        lvalve_length=0.636);
}

exa_clr();
