#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#include "utils.h"

int main()
{
        trexio_exit_code rc;
        double energy;
        int32_t n_up;
        int32_t mo_num;
        int64_t n_integrals;

//--------------------------------------------------------------------------------//
//                                 OPENING FILE                                   //
//--------------------------------------------------------------------------------//

        trexio_t* file = trexio_open("data/h2o.h5",'r',TREXIO_AUTO, &rc);
        if(rc !=TREXIO_SUCCESS)
        {
                printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
                exit(1);
        }


//--------------------------------------------------------------------------------//
//                        READING NUCLEAR REPULSION ENERGY                        //
//--------------------------------------------------------------------------------//

        rc = trexio_read_nucleus_repulsion(file, &energy);

        if (rc != TREXIO_SUCCESS)
        {
                printf("TREXIO Error reading nuclear repulsion energy:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
        }

//--------------------------------------------------------------------------------//
//                        OBTAINING THE NUMBER OF OCCUPIED ORBITALS               //
//--------------------------------------------------------------------------------//

        rc = trexio_read_electron_up_num(file, &n_up);
        if (rc != TREXIO_SUCCESS)
        {
                printf("TREXIO Error reading number of up spin electrons:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
        }

//--------------------------------------------------------------------------------//
//                          READING THE NUMBER OF MOLECULAR ORBITALS              //
//--------------------------------------------------------------------------------//
        rc = trexio_read_mo_num(file, &mo_num);
        if (rc != TREXIO_SUCCESS)
        {
                printf("TREXIO Error reading number of molecular orbitals:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
        }

//--------------------------------------------------------------------------------//
//                          READING ONE-ELECTRON INTEGRALS                        //
//--------------------------------------------------------------------------------//

        double* data = malloc(mo_num * mo_num * sizeof(double));
        if (data == NULL)
        {
                printf("Memory allocation failed for one-electron integrals!\n");
                trexio_close(file);
                return -1;
        }

        rc = trexio_read_mo_1e_int_core_hamiltonian(file, data);
        if (rc != TREXIO_SUCCESS)
        {
                printf("TREXIO Error reading one-electron integrals:\n%s\n", trexio_string_of_error(rc));
                free(data);
                trexio_close(file);
                exit(1);
        }

//--------------------------------------------------------------------------------//
//                          READING TWO-ELECTRON INTEGRALS SIZE                   //
//--------------------------------------------------------------------------------//
        rc = trexio_read_mo_2e_int_eri_size(file, &n_integrals);
        if (rc != TREXIO_SUCCESS)
        {
                printf("TREXIO Error reading two-electron integrals size:\n%s\n",trexio_string_of_error(rc));
                free(data);
                trexio_close(file);
                exit(1);
        }

//--------------------------------------------------------------------------------//
//                          READING TWO-ELECTRON INTEGRALS                        //
//--------------------------------------------------------------------------------//

        int32_t* const index = malloc(4 * n_integrals * sizeof(int32_t));
        if (index == NULL) {
                printf("Memory allocation failed for two-electron integral indices!\n");
                free(data);
                trexio_close(file);
                return -1;
        }

        double* const value = malloc(n_integrals * sizeof(double));
        if (value == NULL) {
                printf("Memory allocation failed for two-electron integral values!\n");
                free(data);
                free(index);
                trexio_close(file);
                return -1;
        }
	
	int64_t buffer_size = n_integrals;

        rc = trexio_read_mo_2e_int_eri(file, 0, &buffer_size, index, value);
        if (rc != TREXIO_SUCCESS)
        {
                printf("TREXIO Error reading two-electron integrals:\n%s\n",trexio_string_of_error(rc));
                free(data);
                free(index);
                free(value);
                trexio_close(file);
                exit(1);
        }

//--------------------------------------------------------------------------------//
//                         READING THE INTEGRAL VALUES                            //
//--------------------------------------------------------------------------------//

	int size_i, size_j, size_k, size_l;

	// Allocation of 4D array with a size of mo_num for each dimension
	double ****integral_array = malloc_4d(mo_num, mo_num, mo_num, mo_num);


	for (int n = 0; n < n_integrals; n++) 
	{
    		int i = index[4 * n + 0];
    		int j = index[4 * n + 1];
    		int k = index[4 * n + 2];
    		int l = index[4 * n + 3];

    		// The integral values considering the 8-fold symmetry
    		integral_array[i][j][k][l] = value[n];
    		integral_array[k][l][i][j] = value[n];
    		integral_array[k][j][i][l] = value[n];
    		integral_array[i][l][k][j] = value[n];
    		integral_array[j][i][l][k] = value[n];
    		integral_array[l][k][j][i] = value[n];
    		integral_array[j][k][l][i] = value[n];
    		integral_array[l][i][j][k] = value[n];
	}

	// Free the allocated memory when done
	free_4d(integral_array, mo_num, mo_num, mo_num);


//--------------------------------------------------------------------------------//
//                         READING THE ORBITAL ENERGIES                           //
//--------------------------------------------------------------------------------//
        
	rc = trexio_read_mo_energy(file, mo_energy);
	if (rc != TREXIO_SUCCESS) 
	{
                printf("TREXIO Error reading molecular orbital energies:\n%s\n",trexio_string_of_error(rc));
                exit(1);
        }

        double mp2_energy = calculate_mp2_energy(mo_energy, index, value, n_up, mo_num, n_integrals);
        printf("MP2 Energy Correction: %.8f Hartree\n", mp2_energy);


//--------------------------------------------------------------------------------//
//                                 CLOSING THE FILE                               //
//--------------------------------------------------------------------------------//
        
        rc = trexio_close(file);
        if (rc != TREXIO_SUCCESS)
        {
                printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
                exit(1);
        }
        file = NULL;

        return 0;
}


//--------------------------------------------------------------------------------//
//                          CALCULATING THE HARTREE FOCK AND MP2 ENERGY           //
//--------------------------------------------------------------------------------//

        double hf_energy = calculate_hartree_fock_energy(energy, data, index, value, n_up, mo_num, n_integrals);
        printf("Hartree-Fock Energy: %.8f Hartree\n", hf_energy);

        double* mo_energy = malloc(mo_num * sizeof(double));
        if (mo_energy == NULL) {
                printf("Memory allocation failed for orbital energies!\n");
                free(data);
                free(index);
                free(value);
                trexio_close(file);
                return -1;
        }
	


