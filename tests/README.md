# Test cases for MFIX-Exa

## Overview of fluid (FLD) cases

***_FLD01:_ 2D Channel (Poiseuille) flow in a periodic domain [1]
***_FLD02:_ 2D Couette flow in a channel [1]
***_FLD03:_ 2D Channel (Poiseuille) flow with specified pressure boundaries  [1]
***_FLD04:_ Tracer transport with and without particles
***_FLD08:_ Flow through a rotated pipe

## Overview of DEM and Fluid/DEM coupled cases

***_DEM01:_ Freely falling particle w/wall collision** [1]

A particle free-falls in a vacuum and collides with a wall. The position and velocity during the three stages (free-fall, collision, rebound) are reported. This case serves to very the linear spring-dashpot collision model as well as the accuracy of the time-stepping methods.


***_DEM02:_ Bouncing particle height** [1]

A particle free-falls in a vacuum and collies with a wall several times. The resulting height the particle reaches post-collision is reported. This case provides a comparison between the linear spring-dashpot model and the hard-sphere model where collisions are instantaneous.


***_DEM03:_ Two stacked, compressed particles** [1]

Two particles are compressed between two walls, and their center positions are reported. This case serves to verify the linear spring-dashpot collision model through analysis of a multi-particle, enduring collision.

***_DEM04:_ Ball slipping on a rough surface** [1]

A spherical particle is rolled on a rough surface. The angular and translational velocities are reported. This case serves to verify the soft-spring collision model through the analysis of the rolling friction model.

***_DEM05:_ Oblique particle collision** [1]

Ninety-three cases of oblique particle collision are tested simultaneously; the particles collide
only with the hard surface, not each other.  Rather than having the particles fall vertically and
hit an angled surface, the wall is kept level (flat) and the particles are given initial trajectories,
resulting in oblique collisions. The particles are initially positioned close to the wall and
gravity is suppressed.  The post-collision angle and angular velocity are reported for each particle.
This case serves to verify the normal and tangential components of both the linear spring-dashpot model.

***_DEM06:_ Single particle terminal velocity** [1]

A particle, initially at rest, is released in a uniformly flow fluid. The particle's velocity increases until it reaches _terminal velocity_ where the gravitational force balances the fluid-particle drag force. The calculated terminal velocity is compared to the theoretical value.


***_DEM07:_ Homogeneous Cooling System** [2]

Particles with no net flow (zero mean velocities) and specified initial granular temperature is allowed to cool over time. The granular temperature of the system is reported over time.


## References
1. J. Musser and A. Choudhary, "MFIX Documentation Volume 3: Verification and Validation Manual," from URL http://mfix.netl.doe.gov

2. Peiyuan Liu, Timothy Brown, William D. Fullmer, Thomas Hauser, and Christine Hrenya, "A Comprehensive Benchmark Suite for Simulations of Particle Laden Flows Using the Discrete Element Method with Performance Profiles from the Multiphase Flow with Interface eXchanges (MFiX) Code," Technical Report NREL/TP-2C00-65637, January 2016.
