#include <stdio.h>
#include <trexio.h>

int main() {
	// opening a trexio file
	trexio_exit_code rc;
	trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
	if (rc != TREXIO_SUCCESS) {
	printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	exit(1);
	}

	//Reading the nuclear repulsion energy
	trexio_exit_code trexio_read_nucleus_repulsion(trexio_t* const trexio_file,double* const energy);

	trexio_exit_code rc; // TREXIO return code
	double energy; // Variable where the energy is read
	rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
	// Check the return code to be sure reading was OK
	if (rc != TREXIO_SUCCESS) {
	printf("TREXIO Error reading nuclear repulsion energy:\n%s\n",
	trexio_string_of_error(rc));
	exit(1);
	}

	//Obtaining the number of occupied orbitals 
	trexio_exit_code trexio_read_electron_up_num(trexio_t* const trexio_file,int32_t* const n_up);

	//Reading one-electron integrals
	trexio_exit_code trexio_read_mo_num(trexio_t* const trexio_file,int32_t* const mo_num);

	trexio_exit_code trexio_read_mo_1e_int_core_hamiltonian(trexio_t* const trexio_file,double* const data);

	//Reading two-electron integrals
	trexio_exit_code trexio_read_mo_2e_int_eri_size(trexio_t* const trexio_file,int64_t* const n_integrals);
		//Allocating memory
	int32_t* const index = malloc(4 * n_integrals * sizeof(int32_t));
	if (index == NULL) {
	fprintf(stderr, "Malloc failed for index");
	exit(1);
	}
	double* const value = malloc(n_integrals * sizeof(double));
	if (value == NULL) {
	fprintf(stderr, "Malloc failed for value");
	exit(1);
	}
		//Reading the integrals
	trexio_exit_code trexio_read_mo_2e_int_eri(trexio_t* const file,const int64_t offset_file,int64_t* const buffer_size,int32_t* const index,double* const value);

		//accesing the n-th integral
	int i = index[4*n+0];
	int j = index[4*n+1];
	int k = index[4*n+2];
	int l = index[4*n+3];
	double integral = value[n];

	//Compute the Hartree-Fock energy
	
	//Compute the MP2 energy
	trexio_exit_code trexio_read_mo_energy(trexio_t* const file,double* const mo_energy);

	//closing a trexio file
	rc = trexio_close(trexio_file);
	if (rc != TREXIO_SUCCESS) {
	printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	exit(1);
	}
	trexio_file = NULL;




