#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ---------------------------------------------------------------------------------------------//
//                      TO ALLOCATE AND FREE THE MEMORY FOR 2-D ARRAY AND MASS                  //
// ---------------------------------------------------------------------------------------------//

double** malloc_2d(size_t m, size_t n)
{
    // Allocate an array of double pointers with size m
    double** a = malloc(m * sizeof(double*));
    if (a == NULL)
    {
        return NULL;
    }

    // Allocate a contiguous block of memory for the 2D array elements
    a[0] = malloc(n * m * sizeof(double));
    if (a[0] == NULL)
    {
        free(a);
        return NULL;
    }

    // Set the pointers in the array of double pointers to point to the correct locations
    for (size_t i = 1; i < m; i++)
    {
        a[i] = a[i - 1] + n;
    }

    return a;
}

void free_2d(double** a)
{
    free(a[0]);
    a[0] = NULL;
    free(a);
}



// ---------------------------------------------------------------------------------------------//
// 			TO READ THE NUMBER OF ATOMS FROM THE INP.TXT FILE                       // 
// ---------------------------------------------------------------------------------------------//

size_t read_Natoms(FILE* input_file)
{
    size_t Natoms;

    // Check if the input file is valid
    if (input_file == NULL)
    {
        printf("Error: Input file is NULL\n");
        return 0;
    }

    // Read the number of atoms from the first line of the file
    if (fscanf(input_file, "%zu", &Natoms) != 1)
    {
        printf("Error: Failed to read the number of atoms\n");
        return 0;
    }

    return Natoms;
}



// ---------------------------------------------------------------------------------------------//
//               TO READ THE COORDINATES AND MASS OF ATOMS FROM THE INP.TXT FILE                //
// ---------------------------------------------------------------------------------------------//

int read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass)
{
    // Check if the input file is valid
    if (input_file == NULL)
    {
        printf("Error: Input file is NULL\n");
        return 0;
    }

    // Read the atomic coordinates and masses for each atom
    for (size_t i = 0; i < Natoms; i++)
    {
        // Read x, y, z coordinates and mass
        if (fscanf(input_file, "%lf %lf %lf %lf", &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]) != 4) {
            printf("Error: Failed to read data for atom %zu\n", i + 1);
            return 0;
        }
    }

    return 1;
}

// ---------------------------------------------------------------------------------------------//
//               		TO CALCULATE THE DISTANCE BETWEEN ATOMS                		//
// ---------------------------------------------------------------------------------------------/

void compute_distances(size_t Natoms, double** coord, double** distance) 
{
    for (size_t i = 0; i < Natoms; i++) 
    {
        for (size_t j = 0; j < Natoms; j++) 
	{
            if (i == j) 
	    {
                distance[i][j] = 0.0;  // Distance from atom to itself
            } 
	    else 
	    {
                double dx = coord[i][0] - coord[j][0];
                double dy = coord[i][1] - coord[j][1];
                double dz = coord[i][2] - coord[j][2];
                
                // Euclidean distance
                distance[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
            }
        }
    }
}
