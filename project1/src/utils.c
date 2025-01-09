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
double calculate_hartree_fock_energy(double nuclear_repulsion_energy, double* one_e_integrals, 
                                     int32_t* indices, double* two_e_values, int32_t n_occ, 
                                     int32_t mo_num, int64_t n_integrals) {
        double hf_energy = nuclear_repulsion_energy;
        double one_e_sum = 0.0;
        double two_e_sum = 0.0;

        // One-electron term
        for (int i = 0; i < n_occ; i++) {
                one_e_sum += one_e_integrals[i * mo_num + i]; // Diagonal terms only
        }
        hf_energy += 2.0 * one_e_sum;

        // Two-electron term
        for (int64_t n = 0; n < n_integrals; n++) {
                int i = indices[4 * n + 0];
                int j = indices[4 * n + 1];
                int k = indices[4 * n + 2];
                int l = indices[4 * n + 3];
                double value = two_e_values[n];

                if (i < n_occ && j < n_occ && k < n_occ && l < n_occ) {
                        two_e_sum += 2.0 * value;  // Coulomb term
                        if (j != k) {  // Avoid double-counting same indices
                                two_e_sum -= value;    // Exchange term
                        }
                }
        }
        hf_energy += two_e_sum;

        return hf_energy;
}

// Function to calculate MP2 energy correction
double calculate_mp2_energy(double* mo_energy, int32_t* indices, double* two_e_values, int32_t n_occ, int32_t mo_num, int64_t n_integrals) 
{
        double mp2_energy = 0.0;

        for (int i = 0; i < n_occ; i++) {
                for (int j = 0; j < n_occ; j++) {
                        for (int a = n_occ; a < mo_num; a++) {
                                for (int b = n_occ; b < mo_num; b++) {
                                        // Locate the integral ⟨ij|ab⟩
                                        double integral_ijab = 0.0;
                                        double integral_ijba = 0.0;

                                        for (int64_t n = 0; n < n_integrals; n++) {
                                                if (indices[4 * n + 0] == i && indices[4 * n + 1] == j &&
                                                    indices[4 * n + 2] == a && indices[4 * n + 3] == b) {
                                                        integral_ijab = two_e_values[n];
                                                }
                                                if (indices[4 * n + 0] == i && indices[4 * n + 1] == j &&
                                                    indices[4 * n + 2] == b && indices[4 * n + 3] == a) {
                                                        integral_ijba = two_e_values[n];
                                                }
                                        }

                                        // Compute MP2 contribution
                                        double denominator = mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b];
                                        mp2_energy += (integral_ijab * (2.0 * integral_ijab - integral_ijba)) / denominator;
                                }
                        }
                }
        }

        return mp2_energy;
}
