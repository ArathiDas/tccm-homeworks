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
	
    // Read the number of atoms from the first line of the file
    if(fscanf(input_file, "%zu", &Natoms) != 1)
    {
	    return 0;
    }
    
    return Natoms;
}



// ---------------------------------------------------------------------------------------------//
//               TO READ THE COORDINATES AND MASS OF ATOMS FROM THE INP.TXT FILE                //
// ---------------------------------------------------------------------------------------------//

int read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass)
{
    
    // Read the atomic coordinates and masses for each atom
    for (size_t i = 0; i < Natoms; i++)
    {
        // Read x, y, z coordinates and mass
        if(fscanf(input_file, "%lf %lf %lf %lf", &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]) != 4)
	{
		return 0;
	}
    }

    return 1;
}

// ---------------------------------------------------------------------------------------------//
//               		TO CALCULATE THE DISTANCE BETWEEN ATOMS                		//		
// ---------------------------------------------------------------------------------------------//

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
                // The distances are multiplied by 0.1 unit to make the "nm" scale
		double dx = 0.1 * (coord[i][0] - coord[j][0]);
                double dy = 0.1 * (coord[i][1] - coord[j][1]);
                double dz = 0.1 * (coord[i][2] - coord[j][2]);
                
                // Euclidean distance
                distance[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------//
//                              TO CALCULATE THE POTENTIAL ENERGY                               //
// ---------------------------------------------------------------------------------------------//

double potential_energy(double epsilon, double sigma, size_t Natoms, double** distance)
{
	double V_LJ = 0.0;
	double V_total = 0.0;

	for (int i=0; i<Natoms; i++)
	{
		for (int j=i+1; j<Natoms; j++)
		{
			double r = distance[i][j];
			double sigma_over_r_value;  
			double power_12_term_value; 			
			double power_6_term_value;
			
			if (r > 0)
			{
				sigma_over_r_value  = sigma/r;
				power_12_term_value = pow(sigma_over_r_value,12);
				power_6_term_value  = pow(sigma_over_r_value,6);
				
				V_LJ = 4*epsilon*(power_12_term_value - power_6_term_value);
			}

			V_total += V_LJ;
			V_LJ = 0.0;

		}
	}
	return V_total;

}
					

	
