
#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>

// Function to allocate the memory
double** malloc_2d(size_t m, size_t n);

// Functions to free the allocated memory
void free_2d(double** a);

// Function to read the number of atoms from the input file (inp.txt)
size_t read_Natoms(FILE* input_file);

// Functions to read the coordinate and mass of the molecule from the inpiut file (inp.txt)
int read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass);

// Function to compute the internuclear distance between pairs of atoms
void compute_distances(size_t Natoms, double** coord, double** distance);

// Function to compute the potential energy 
double potential_energy(double epsilon, double sigma, size_t Natoms, double** distance);
#endif

