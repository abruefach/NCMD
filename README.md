# Nanocluster Molecular Dynamics (NCMD)
## This project was initially developed for University of California's Computational Nanomechanics Class (CE237) in Spring 2021.

This is an interactive tutorial for demonstrating molecular dynamics of nanocluster systems using Langevin Dynamics in python. We implement methods from Atomic Simulation Environment (ASE) and As Soon As Possible (ASAP). 

The growth process is modelled from stage II to III of the LaMer Curve (1). We make the simlification that the system is in equilibrium, thus allowing us to use simple dynamics to simulate the growth. While not entirely accurate, this model can give us insights on how nuclei form in solution.

## Simulation Details

Potential - All potential energy calculators in ASE or ASAP are supported in this workflow. Our provided example uses the Lennard Jones potential.

Velocity Distribution - This approach currently uses the Maxwell Boltzmann Distribution.

Temperature Control - We use Langevin Thermostat in this example, which is an equilibrium approach at constant NVT. The Langevin approach holds temperature constant by applying a drag force to mitigate large temperature spikes.

## For Users

You can download the files in this repository by using the git command or by downloading the .zip folder. Make sure the dependencies are installed to use.

## Dependencies

- Atomic Simulation Environment (ASE)
- As soon as possible (ASAP)
- Numpy
- Scipy
- Matplotlib

## References
1. Pound, Guy M., and Victor K. La Mer. "Kinetics of crystalline nucleus formation in supercooled liquid Tin1, 2." Journal of the American Chemical Society 74.9 (1952): 2323-2332.

## License
GNU GPLv3

NCMD is distributed under a GPLv3 license. It is free to use, alter, or build on, provided that any work derived from NCMD is also kept free and open.
