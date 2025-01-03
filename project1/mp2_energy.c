#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>

int main()
{
	trexio_exit_code rc;
	double energy;
	int32_t n_up;
	int32_t mo_num;

//--------------------------------------------------------------------------------//
//				   OPENING FILE				          // 
//--------------------------------------------------------------------------------//

	trexio_t* file = trexio_open("h2o.h5",'r',TREXIO_AUTO, &rc);
	if(rc !=TREXIO_SUCCESS)
	{
		printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
		exit(1);
	}


//--------------------------------------------------------------------------------//
//              	  READING NUCLEAR REPULSION ENERGY                        //
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

        // Allocate memory for the data array
        double* data = (double*)malloc(mo_num * mo_num * sizeof(double));
        if (data == NULL)
        {
                printf("Memory allocation failed!\n");
                return -1;
        }

        rc = trexio_read_mo_1e_int_core_hamiltonian(file, data);
        if (rc != TREXIO_SUCCESS)
        {
                printf("TREXIO Error reading number of up spin electrons:\n%s\n",trexio_string_of_error(rc));
                exit(1);
        }

//--------------------------------------------------------------------------------//
//                          READING TWO-ELECTRON INTEGRALS SIZE                   //
//--------------------------------------------------------------------------------//
	rc = trexio_read_mo_2e_int_eri_size(file,&n_integrals);
	if (rc != TREXIO_SUCCESS)
        {
        	printf("TREXIO Error reading two-electron integrals size:\n%s\n",trexio_string_of_error(rc));
                exit(1);
        }

//--------------------------------------------------------------------------------//
//                          ALLOCATING MEMORY FOR TWO-ELECTRON INTEGRALS          //
//--------------------------------------------------------------------------------//
	int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
 	if (index == NULL) {
 		fprintf(stderr, "Malloc failed for index");
 		exit(1);
	}
 	double* value = malloc(n_integrals * sizeof(double));
 		if (value == NULL) {
 		fprintf(stderr, "Malloc failed for value");
 		exit(1);
 	}

//--------------------------------------------------------------------------------//
//                          READING TWO-ELECTRON INTEGRALS                   //
//--------------------------------------------------------------------------------//
	int64_t buffer_size = n_integrals;
	int64_t offset_file = 0;
	rc = trexio_read_mo_2e_int_eri(file, offset_file, &buffer_size, index, value);
	if (rc != TREXIO_SUCCESS)
        {
        	printf("TREXIO Error reading two-electron integrals:\n%s\n",trexio_string_of_error(rc));
                exit(1);
        }
	
//--------------------------------------------------------------------------------//
//                              PRINTING THE INPUT VALUE                          //
//--------------------------------------------------------------------------------//


    printf("Nuclear Repulsion Energy: %.6f Hartree\n", energy);
    printf("Number of up spin electrons: %d \n", n_up);
    printf("Number of molecular orbitals: %d \n", mo_num);

    for(int i=0; i<mo_num; i++)
    {
            for(int j=0; j<mo_num; j++)
            {
                    printf("The one electron integral [%d][%d]=%f \n ",i,j,data[i*mo_num + j]);
            }
    }

	
//--------------------------------------------------------------------------------//
//                                 CLOSING FILE                                   //
//--------------------------------------------------------------------------------//    
	rc = trexio_close(file);
	
	if (rc != TREXIO_SUCCESS) 
	{
		printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	file = NULL;

}




