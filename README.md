# DNA–Microtubule Network Simulation

Flexible spring–bead polymer model ("DNA") in a 2D square box with periodic boundary conditions (PBC). Bonds use FENE springs; non‑bonded interactions switch from a Gaussian soft volume repulsion during early iterations to Lennard–Jones (LJ) afterwards. The terminal bead of each polymer experiences active self‑propulsion with an anisotropic friction tensor. Thermal noise uses Gaussian random numbers from Intel MKL VSL. Parallelism is via OpenMP.

This repository contains the simulation code used for the study:
> "Dynamic assembly of complex hierarchical DNA polymer networks by biomolecular active agents"
> and partially based on the simulation done by Rakesh Das, et.al. (2022) "How enzymatic activity is involved in chromatin organization, eLife 11:e79901."

----
## Features
- **Polymer model**: FENE(Finite Extensible Non-linear Elastic) + Lennard–Jones potentials
- **Self-propulsion**: Applied to the end monomer with defined polarity
- **Periodic boundary conditions** in 2D
- **Example configuration** provided for test runs

##Compilation
-  Load MKL (adjust to your cluster)
module load mkl/2022.1.0

# Compile
gfortran DNA-MT.f90 -o DNA-MT.exe \
  -I"${MKLROOT}/include" \
  -L"${MKLROOT}/lib/intel64" \
  -Wl,--start-group -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group \
  -fopenmp -lpthread -lm -ldl -O3

If you see undefined reference to vsl... errors, double‑check the include path, library path, MKL link order, and that -fopenmp is present.

##Requirements
Compiler: GNU Fortran (gfortran) with OpenMP support.
Intel MKL (tested with mkl/2022.1.0) — VSL headers and libraries are required (mkl_vsl.f90, VSL symbols).
Pthreads, libm, libdl
Module environment (optional): Example uses module load mkl/2022.1.0.
#All input files are in input-files.zip, and are in free format; change the iteration value for shorter runs. Counts must be consistent.
#An example of 8 polymers, consisting of 64 monomers each, is attached in the input-pos.inp for test runs.
