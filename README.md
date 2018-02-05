# Trappist-1 Stellar System Simulation

system parameters source:
http://www.trappist.one/#system

This goal of this project was to create a N-body gravitational motion simulator. The motivation behind this project came
from the discovery of the Trappist-1 exoplanetary system, which is known to contain about seven planets of Earth-like mass
all orbiting within .07 AU of their host star. This system is much more crowded than our own solar system would seem 
perhaps unstable and I wondered how resilient it would be to changes of the any of the planetary masses. Significant uncertainties
still  exist in measurements of the planets' masses, and perhaps and suitable simulation of the system could rule out configurations
that result in instability. 

The create the simulation, the differential equation solver program using the Euler algorithm from Session 6 was adapted so that it 
could solve the second order gravitational acceleration differential equation in multiple dimensions and for multiple objects. It is 
necessary for this numerical calculation of the equations of motion, since each object's acceleration depends on the position of 
every other object, which are also themselves moving and so a solution by simple mathematical manipulation is not possible and numerical
calculations must be performed at regular time intervals (the mesh size of Euler alg.). 

The program was tested at incremental stages of development, first starting with just simulations of the Moon orbiting Earth to confirm
expected motions and then by test cases made by altering masses and distances. Eventually limitations of runtime and data file size were 
identified. It was desirable to simulate over the course of many periods to observe long term trends in the orbital motions and for the 
Trappist-1 system, this meant about 100 days (8e6 seconds) for 10 periods of the furthest planets. Dimension was limited to 2 (planets
approximately orbit in the same plane) and mesh size was taken to be 2 seconds to yield a runtime of about 20 minutes and data files of
333 MB for each body. These parameters could be increased if more time/computational power/data storage is available. Once the data files for each body were generated, a script loop created gnuplot images of the system at incremental times, along with their past trajectory. Finally, each of these images could be loaded into gimp to create an animation of the motion for easier vizualization and interpretation.

The planets were assumed to start in complete alignment with a horizontal velocity given by a circular orbit (which is approximately true) from the equation v = sqrt(GM\star/r). Also the most probable value of the planets masses were used (meaning the value was chosen by neglecting uncertainty) With a mesh size of 2 seconds, the system's inner planets were seen to "drift" in the positive x direction of the course of the 100 days. It is uncertain if this was due to system instability or to global error in the Euler algorithm, or to the choice of initial conditions. However, a second run was attempted, with a larger mesh size of 10 seconds and interestingly, the drift effect was larger, but in the same direction. This suggests that the error from the mesh size plays an important role, but that there is another cause, or else we would expect the error to point in a different random direction. Most likely it is due to the initial tug from all of the planets starting in alignment. The trend of Trappist-B nearing the orbital path of Trappist-C suggests that long term this is not a stable configuration and that in fact all initial conditions are NOT viable. This says the system is in fact sensitive to perturbations and more testing could be done with different initial conditions and masses for each planet to try to determine if there are any more stable configurations and also if over time all simulations will be doomed to instability my accumulating global errors of the solution algorithm. Ideally, the program would test over a wide range of parameters, however this is computationally intensive and a more optimized solution would need to be found. 

Several additonal tests were created with a mesh 10 - firstly the lowest mass estimate of planet B was used, while the highest estimate of C and D were used. Following this, the masses of C and D were upgraded to over twice that of their upper estimate. This was found to increase the effect the planetary drift, but not by an overly significant amount. These tests suggest than much longer simulation time would be needed to observe the full effects of changing mass parameters. Conversely, however, the initial locations of planets C-G were flipped over the y-axis and this was found to completely reverse the direction of planet B's drift. This lends support to the idea that the system's final state is not independent of its initial conditions. This only has limited relevence to a planetary system, as the planets were not one day placed at a given spot with an initial velocity, but instead gradually formed from accretion. It still does suggest  that the given orientation can have a dramatic impact on the time evolution of the system and contributions from small perturbations can add up. 


