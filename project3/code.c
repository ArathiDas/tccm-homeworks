#include <stdio.h>
#include <math.h>

//ALLOCATING AND FREEING 2D ARRAYS

double** malloc_2d(size_t m, size_t n) {
  // Allocate an array of double pointers with size m
  double** a = malloc(m*sizeof(double*));
  if (a == NULL) {
  	return NULL;
  }
  // Allocate a contiguous block of memory for the 2D array elements
  a[0] = malloc(n*m*sizeof(double));
  if (a[0] == NULL) {
  	free(a);
	return NULL;
  }
  // Set the pointers in the array of double pointers
  // to point to the correct locations
  for (size_t i=1 ; i<m ; i++) {
	 a[i] = a[i-1]+n;
  }
  return a;
}

void free_2d(double** a) {
	free(a[0]);
        a[0] = NULL;
	free(a);
}


//DESCRIBING THE ATOMS

size_t read_Natoms(FILE* input_file);

void read_molecule(FILE* input_file,
	size_t Natoms,
	double** coord,
	double* mass);

void compute_distances(size_t Natoms,
	double** coord,
	double** distance);

// LENNARD - JONES POTENTIAL 
double V(double epsilon,
	double sigma,
	size_t Natoms,
	double** distance);

// COMPUTING THE TOTAL ENERGY
double T(size_t Natoms,
	double** velocity,
	double* mass);

// COMPUTING THE ACCELERATION
void compute_acc(size_t Natoms,
	double**coord,
	double* mass,
        double**distance,
        double**acceleration);

// IMPLEMENTING THE MOLECULAR DYNAMICS 
