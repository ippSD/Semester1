# orbits

Orbital Dynamics library. Includes:

* N-Bodies & Kepler functions according to the Cauchy Problem style (orbit_functions.f90).
* Object Oriented classes (orbit_objects.f90).
* Lagrange Points of two-bodies systems (orbit_lagrange_points.f90)
* Subroutines for setting pointers to the state vector for easily accessing
the required data without messing with state vector indexes (orbit_pointers.f90).
* Checks for asserting correct inputs (orbit_check.f90).
* Subroutines for exporting orbit data to data files (orbit_exports.f90).
* Quick orbit plots (orbit_plots.f90).