#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>

int main()
{
	trexio_exit_code rc;
	double energy;
	int32_t n_up;

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
//                        READING NUCLEAR REPULSION ENERGY                        //
//--------------------------------------------------------------------------------//

	rc = trexio_read_electron_up_num(file, &n_up);
	if (rc != TREXIO_SUCCESS)
        {
                printf("TREXIO Error reading Number of up spin electrons:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
        }


//--------------------------------------------------------------------------------//
//                              PRINTING THE INPUT VALUE	                  //                   
//--------------------------------------------------------------------------------//


    printf("Nuclear Repulsion Energy: %.6f Hartree\n", energy);
     printf("Number of up spin electrons: %d \n", n_up);




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
                printf("TREXIO Error reading Number of up spin electrons:\n%s\n",trexio_string_of_error(rc));
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


