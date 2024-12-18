#include<stdio.h>
#include<stdlib.h>
#include "utils.h"
#include "error.h"
// ---------------------------------------------------------------------------------------------//
//                      		THE MAIN PROGRAM    					//
// ---------------------------------------------------------------------------------------------//


int main() 
{
    FILE* input_file = fopen("inp.txt", "r"); 			// Open the file in read mode
    error_file_open(input_file);

    // Read the number of atoms in the molecule
    size_t Natoms = read_Natoms(input_file);
    if (Natoms == 0) 
    {
        error_read_atoms();
    }

    // --------------------- Allocation and initialization of the 2-D array -------------------------
    double** coord = malloc_2d(Natoms, 3);
    if (coord == NULL) 
    {
       error_memory_allocation(coord);
    }

    double* mass = malloc(Natoms * sizeof(double)); 		// Mass array for atoms
    if (mass == NULL) 
    {
        error_memory_allocation(mass);
        free_2d(coord); 					// Free previously allocated memory
          
    }

    double** distance = malloc_2d(Natoms,Natoms);
    if (distance == NULL)
    {
       error_memory_allocation(coord);
    }

    // Read the molecule data from the file
    if (read_molecule(input_file, Natoms, coord, mass) != 1) 
    {
        free_2d(coord);    					// Free allocated memory
        free(mass);
        fclose(input_file); 					// Close the file
        return 1;
    }


    // ---------------- Calculating the internuclear distance using the coordinates -----------------
    compute_distances(Natoms,coord,distance);



    // ----------------------- Print the data to verify it was read correctly ------------------------
    for (size_t i = 0; i < Natoms; i++) 
    {
        printf("Atom %zu: x = %lf, y = %lf, z = %lf, mass = %lf\n",
               i + 1, coord[i][0], coord[i][1], coord[i][2], mass[i]);
    }
    
   
    printf("Internuclear Distances:\n");
    for (size_t i = 0; i < Natoms; i++) {
        for (size_t j = 0; j < Natoms; j++) {
            printf("d(%zu,%zu) = %lf ", i + 1, j + 1, distance[i][j]);
        }
        printf("\n");
    }

    free_2d(coord);     					// To free the memory 2-D array "coord"
    free(mass);         					// To free the memory array "mass"
    free_2d(distance);
    fclose(input_file); 					// Close the file

    return 0;
}

