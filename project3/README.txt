# MOLECULAR DYNAMICS SIMULATION IN C

This project implements a Molecular Dynamics (MD) simulation in C. The code simulates the motion of the atoms using Lennard-Jones Potential, computes kinetic and potential energy and outputs atomic trajectories for visualization. 

## Features 
- Reads atomic data from an input file 
- Computes internuclear distance and Lennard-Jones Potential
- Implements the Verlet algorithm for time integration
- Outputs atomic trajectories in XYZ format for visualization with tools like Molden

## Directory structure
- LICENSE: contains information about the terms and conditions for the usage and redistribution of the code
- AUTHORS: Contains the names of individuals who contributed to the development of the project
- INSTALL.md: Provides instructions on how to compile and run the program
- Makefile: handles the compilation process for the source files
- data/: contains the input files for the program for methane, water and benzene
- src/: Source code for the simulation 
     - dynamics.c: implements the core dynamics
     - utils.c: contains utility functions for memory allocation, reading inputs, defining functions, etc.
     - utils.h: header file for utility functions declarations
     - error.c: defines error-handling functions for the program
     - error.h: header file for error-handling functions
