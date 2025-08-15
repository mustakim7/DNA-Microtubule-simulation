# DNA–Microtubule Network Simulation

Flexible spring–bead polymer model ("DNA") in a 2D square box with periodic boundary conditions (PBC). Bonds use FENE springs; non‑bonded interactions switch from a Gaussian soft volume repulsion during early iterations to Lennard–Jones (LJ) afterwards. The terminal bead of each polymer experiences active self‑propulsion with an anisotropic friction tensor. Thermal noise uses Gaussian random numbers from Intel MKL VSL. Parallelism is via OpenMP.

This repository contains the simulation code used for the study:
> "Dynamic assembly of complex hierarchical DNA polymer networks by biomolecular active agents"

----
## Features
- **Polymer model**: FENE(Finite Extensible Non-linear Elastic) + Lennard–Jones potentials
- **Self-propulsion**: Applied to end monomer with defined polarity
- **Periodic boundary conditions** in 2D
- **Example configuration** provided for test runs
