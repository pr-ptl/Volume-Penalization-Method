Penalization methods, also known as the Volume Penalization Method (VPM), a technique derived from the broader class of immersed boundary methods (IBM). They offer an effective framework for handling problems involving multiphase flows, internal solid regions, moving boundaries, or complex geometries in fluid and transport simulations.

Rather than explicitly meshing or conforming to the geometry of obstacles or interfaces, VPM embeds both the fluid and solid domains into a larger computational domain, often using a regular grid. The method introduces a penalization term into governing equations (e.g., advection–diffusion or Navier–Stokes), which acts to enforce Dirichlet or Neumann boundary conditions within obstacle regions by damping the solution toward a prescribed target. This is achieved via a mask function that distinguishes fluid from solid regions, and a penalization parameter that controls the enforcement strength.

The functionality of the methods includes
- Simple implementation on structured or unstructured grids
- Efficient when using fast solvers and avoiding complex mesh generation
- Flexible for handling irregular or time-dependent geometries
 - Applicable to a wide variety of physics problems, including EHD, MHD, and diffusion

About This Repository

This repository contains MATLAB scripts demonstrating the Volume Penalization Method for solving the 1D diffusion equation on a regularly spaced domain [−2,2]. The setup imposes Neumann boundary conditions at 
x=-1 and x=1, symmetry conditions at x=2 and x=−2, following the framework described in Kadoch et al. (2012).

The implementation uses:

A first-order Finite Volume Method (FVM) for spatial discretization
An exponential time marching scheme.

A focus on simplicity and clarity, making it ideal for educational use and validation purposes

REFERENCES

[1] Angot, P., Bruneau, C.-H., & Fabrie, P. (1999). A penalisation method to take into account obstacles in viscous flows. Numerische Mathematik, 81(4), 497–520. https://doi.org/10.1007/s002110050357

[2] Kadoch, B., Kolomenskiy, D., Angot, P., & Schneider, K. (2012). A volume penalization method for incompressible flows and scalar advection–diffusion with moving obstacles. Journal of Computational Physics, 231(12), 4365–4383. https://doi.org/10.1016/j.jcp.2012.01.036
