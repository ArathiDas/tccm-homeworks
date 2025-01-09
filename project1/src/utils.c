#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>

// Function to allocate a 4D array to acces the nth integral
double ****malloc_4d(int size_i, int size_j, int size_k, int size_l) 
{    
    // For the first dimension
    double ****array = (double ****)malloc(size_i * sizeof(double ***));
    if (array == NULL) exit(EXIT_FAILURE);

    // For the second dimension
    for (int i = 0; i < size_i; i++) 
    {
        array[i] = (double ***)malloc(size_j * sizeof(double **));
        if (array[i] == NULL) exit(EXIT_FAILURE);

        // For the third dimension
        for (int j = 0; j < size_j; j++) 
	{
            array[i][j] = (double **)malloc(size_k * sizeof(double *));
            if (array[i][j] == NULL) exit(EXIT_FAILURE);

            // For the fourth dimension
            for (int k = 0; k < size_k; k++) 
	    {
                array[i][j][k] = (double *)malloc(size_l * sizeof(double));
                if (array[i][j][k] == NULL) exit(EXIT_FAILURE);
            }
        }
    }

    return array;
}

// Function to free the 4D array
void free_4d(double ****array, int size_i, int size_j, int size_k) 
{
    for (int i = 0; i < size_i; i++) 
    {
        for (int j = 0; j < size_j; j++) 
	{
            for (int k = 0; k < size_k; k++) 
	    {
                free(array[i][j][k]); 
            }
            free(array[i][j]); 
        }
        free(array[i]); 
    }
    free(array); 
}

// Function to calculate Hartree-Fock energy
double HF_energy(double energy, double *data, double ****integral_array, int mo_num, int n_up) 
{
    double one_e_term = 0.0;
    double two_e_term = 0.0;

    // ONE-ELECTRON TERM
    for (int i = 0; i < n_up; i++) 
    {
        one_e_term += data[i * mo_num + i];               // Sum diagonal terms for occupied orbitals
    }
   
    // TWO-ELECTRON TERM
    for (int i = 0; i < n_up; i++) 
    {
        for (int j = 0; j < n_up; j++) 
	{
            two_e_term += 2 * integral_array[i][j][i][j] - integral_array[i][j][j][i];
        }
    }
 
    // HF FINAL ENERGY
    double hf_energy = energy + 2.0 * one_e_term + two_e_term;

    return hf_energy;
}

	
	
// Function to calculate MP2 energy correction
double calculate_MP2_energy(double ****integral_array, double *mo_energy, int n_up,int mo_num) 
{
    double energy_mp2 = 0.0;

    // Loops through occupied orbitals i, j and virtual orbitals a, b
    for (int i = 0; i < n_up; i++) 
    {
        for (int j = 0; j < n_up; j++) 
	{
            for (int a = n_up; a < mo_num; a++) 
	    {
                for (int b = n_up; b < mo_num; b++) 
		{
                    // Compute MP2 energy correction
                    double numerator = integral_array[i][j][a][b] * (2.0 * integral_array[i][j][a][b] - integral_array[i][j][b][a]);
                    double denominator = mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b];
                    energy_mp2 += numerator / denominator;
                }
            }
        }
    }

    return energy_mp2;
}
