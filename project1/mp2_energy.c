#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>

int main()
{
	trexio_exit_code rc;
	double energy;

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
//                              PRINTING THE INPUT VALUE	                  //                   //--------------------------------------------------------------------------------//


    printf("Nuclear Repulsion Energy: %.6f Hartree\n", energy);





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


