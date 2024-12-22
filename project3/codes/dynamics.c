#include<stdio.h>
#include<stdlib.h>
#include "utils.h"
#include "error.h"

// ---------------------------------------------------------------------------------------------//
//                      		THE MAIN PROGRAM    					//
// ---------------------------------------------------------------------------------------------//


int main() 
{
    FILE* input_file = fopen("data/CH4_cluster_1.txt", "r"); 		// Open the file in read mode
    if (input_file == NULL)
    {
	    error_file_open(input_file);			// Calling the error display function
    }

     //-------------------------Read the number of atoms in the molecule-------------------------------//
    
    size_t Natoms = read_Natoms(input_file);			// Calling the read_Natom function to read the number of atom from the input file
    if (Natoms == 0) 
    {
        error_read_atoms();					// Calling the error display function
    }

    // --------------------- Allocation and initialization of the 2-D array ---------------------------//
   
    double** coord = malloc_2d(Natoms, 3);			// Allocation of the 2-D array of size (number of atom as rows and 3 as coordinates) using malloc_2d function 
    if (coord == NULL) 
    {
       error_memory_allocation("coord");			// Calling the error display function
    }

    double* mass = malloc(Natoms * sizeof(double)); 		// Allocation of the Mass array for atoms as size of number of atoms
    if (mass == NULL) 
    {
        error_memory_allocation("mass");			// Calling the error display function
        free_2d(coord); 					// Free function to free up the allocated memory
          
    }

    double** distance = malloc_2d(Natoms,Natoms);		// Allocation of an array "distance" to store the distances, size number of atom by number of atom 2-D array
    if (distance == NULL)
    {
       error_memory_allocation("coord");			// Calling the error display function
    }

    //------------------------Read the molecule data from the file-------------------------------------//
   
    if (read_molecule(input_file, Natoms, coord, mass) != 1) 	// Function to read the coordinates from the input file, this function read only coordinates
    {
        free_2d(coord);    					// Free allocated memory
        free(mass);						// Free allocated mass
        fclose(input_file); 					// Close the file
        return 1;
    }


    // ---------------- Calculating the internuclear distance using the coordinates -------------------//
   
    compute_distances(Natoms, coord, distance);			// Function compute the distances

    // ---------------- Calculating the total potential energy of the sytem------- --------------------//
    
    double epsilon = 0.0661;
    double sigma = 0.3345;
    double V_total;
    V_total = potential_energy(epsilon, sigma, Natoms, distance);


    // ----------------------- Print the data to verify it was read correctly -----------------------//

    for (size_t i = 0; i < Natoms; i++) 
    {
        printf("Atom %zu: x = %lf, y = %lf, z = %lf, mass = %lf\n", i + 1, coord[i][0], coord[i][1], coord[i][2], mass[i]);
    } 
    
   
    printf("Internuclear Distances:\n");
    for (size_t i = 0; i < Natoms; i++) {
        for (size_t j = 0; j < Natoms; j++) {
            printf("d(%zu,%zu) = %lf ", i + 1, j + 1, distance[i][j]);
        }
        printf("\n");
    }
    
    printf("Total Potential Energy: %.6f\n", V_total);

    free_2d(coord);     					// To free the memory 2-D array "coord"
    free(mass);         					// To free the memory array "mass"
    free_2d(distance);
    
    fclose(input_file); 					// Close the file

    return 0;
}


