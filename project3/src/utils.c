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

double V(double epsilon, double sigma, size_t Natoms, double** distance)
{
	double V_LJ = 0.0; 					// Variable to store the Lennard-Jones potential for a pair of atoms
	double V_total = 0.0; 					// Variable to accumulate the total potential energy for all pairs

	// Loop over all unique pairs of atoms
	for (int i=0; i<Natoms; i++)
	{
		for (int j=i+1; j<Natoms; j++) 			// Only consider j > i to avoid double counting
		{
			double r = distance[i][j]; 		// Distance between atom i and atom j
			double sigma_over_r_value;  		// Ratio of sigma to r
			double power_12_term_value; 		// Term for the 1/r^12 contribution
			double power_6_term_value;  		// Term for the 1/r^6 contribution
			
			if (r > 0) // Ensure distance is positive to avoid division by zero
			{
				sigma_over_r_value  = sigma/r;  			// Calculate sigma/r
				power_12_term_value = pow(sigma_over_r_value,12); 	// (sigma/r)^12 term
				power_6_term_value  = pow(sigma_over_r_value,6);  	// (sigma/r)^6 term
				
				// Calculate Lennard-Jones potential for this pair
				V_LJ = 4*epsilon*(power_12_term_value - power_6_term_value);
			}

			V_total += V_LJ;			 // Accumulate the potential energy for this pair
			V_LJ = 0.0; 				 // Reset pair potential for next iteration
		}
	}
	return V_total; 					 // Return the total potential energy for all pairs
}

// ---------------------------------------------------------------------------------------------//
//                              TO CALCULATE THE KINETIC ENERGY                                 //
// ---------------------------------------------------------------------------------------------//

double T(size_t Natoms, double** velocity, double* mass)                     // function to calculate the total kinetic energy  
{
	double kinetic_energy = 0.0;
	double total_kinetic_energy = 0.0;

	for (size_t i = 0; i < Natoms; i++)                                  // iterating over each atom
	{
		double vx = velocity[i][0];                                  // extracting the x, y and z components from the array
		double vy = velocity[i][1];
		double vz = velocity[i][2];

		kinetic_energy = 0.5 *mass[i]*(vx*vx + vy*vy + vz*vz);      // calculating the kinetic energy
		total_kinetic_energy += kinetic_energy;
	}
return total_kinetic_energy;
}

// ---------------------------------------------------------------------------------------------//
//                              TO CALCULATE THE TOTAL ENERGY                                   //
// ---------------------------------------------------------------------------------------------//

double E(double V, double T)
{
	return T+V;
}

// ---------------------------------------------------------------------------------------------//
//                              COMPUTING THE ACCELERATION                                      //
// ---------------------------------------------------------------------------------------------//

// Function to compute accelerations
void compute_acc(size_t Natoms, double** coord, double* mass, double** distance, double** acceleration) {
    // Initialize acceleration to 0
    for (size_t i = 0; i < Natoms; i++) {
        for (size_t j = 0; j < 3; j++) {
            acceleration[i][j] = 0.0;
        }
    }

    // Compute forces and update accelerations
    for (size_t i = 0; i < Natoms; i++) {
        for (size_t j = 0; j < Natoms; j++) {
            if (i != j) {
                double r = distance[i][j];
                if (r > 0) {
                    double sigma_over_r = sigma / r;
                    double power_12 = pow(sigma_over_r, 12);
                    double power_6 = pow(sigma_over_r, 6);
                    double force = (24.0 * epsilon / r) * (2 * power_12 - power_6);

                    acceleration[i][0] += force * (coord[i][0] - coord[j][0]) / (r*mass[i]);
                    acceleration[i][1] += force * (coord[i][1] - coord[j][1]) / (r*mass[i]);
                    acceleration[i][2] += force * (coord[i][2] - coord[j][2]) / (r*mass[i]);
                }
            }
        }
    }
}


// ---------------------------------------------------------------------------------------------//
//                              IMPLEMENTING VERLET ALGORITHM                                    //
// ---------------------------------------------------------------------------------------------//


// Function to implement the Verlet algorithm for molecular dynamics
void verlet_algorithm(size_t Natoms, double dt, size_t steps, 
                      double** coord, double** velocity, 
                      double** acceleration, double** distance, 
                      double* mass, double epsilon, double sigma, 
                      const char* output_file) {
    // Validate inputs
    if (Natoms == 0 || dt <= 0 || steps == 0) {
        printf("Error: Invalid input parameters. Check Natoms, dt, and steps.\n");
        return;
    }

    FILE* outfile = fopen(output_file, "w");
    if (!outfile) {
        printf("Error: Unable to open output file: %s\n", output_file);
        return;
    }

    for (size_t step = 0; step < steps; step++) {
        // Step 1: Update positions
        for (size_t i = 0; i < Natoms; i++) {
            for (size_t j = 0; j < 3; j++) {
                coord[i][j] += velocity[i][j] * dt + 0.5 * acceleration[i][j] * dt * dt;
            }
        }

        // Step 2: Compute new distances
        compute_distances(Natoms, coord, distance);

        // Step 3: Compute new accelerations
        double** new_acceleration = malloc_2d(Natoms, 3); // Temporary storage for new accelerations
        compute_acc(Natoms, coord, mass, distance, new_acceleration);

        // Step 4: Update velocities using old and new accelerations
        for (size_t i = 0; i < Natoms; i++) {
            for (size_t j = 0; j < 3; j++) {
                velocity[i][j] += 0.5 * (acceleration[i][j] + new_acceleration[i][j]) * dt;
                acceleration[i][j] = new_acceleration[i][j]; // Update old acceleration to new
            }
        }

        // Free temporary new_acceleration array
        free_2d(new_acceleration);

        // Step 5: Compute energies
        double potential_energy = V(epsilon, sigma, Natoms, distance);
        double kinetic_energy = T(Natoms, velocity, mass);
        double total_energy = E(potential_energy, kinetic_energy);

        // Write trajectory every 10 steps
        if (step % 10 == 0) {
            fprintf(outfile, "%zu\n", Natoms);  // Number of atoms
            fprintf(outfile, "Step %zu | Kinetic Energy: %.6f | Potential Energy: %.6f | Total Energy: %.6f\n",
                    step, kinetic_energy, potential_energy, total_energy);

            for (size_t i = 0; i < Natoms; i++) {
                fprintf(outfile, "Ar %.6f %.6f %.6f\n", coord[i][0], coord[i][1], coord[i][2]);
            }
        }

        // Energy conservation check and debug log
        printf("Step %zu | Total Energy: %.6f | Kinetic: %.6f | Potential: %.6f\n", 
               step, total_energy, kinetic_energy, potential_energy);
    }

    fclose(outfile);
}




	
