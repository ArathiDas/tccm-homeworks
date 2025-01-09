#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>

//Functon to allocate the memorey 
double ****malloc_4d(int size_i, int size_j, int size_k, int size_l);

//Function to free the memory
void free_4d(double ****array, int size_i, int size_j, int size_k);

//Function reading nuclear repulsion. 
trexio_exit_code trexio_read_nucleus_repulsion(trexio_t* const trexio_file, double* const energy);

//Function reading occupied molecular orbitals
trexio_exit_code trexio_read_electron_up_num(trexio_t* const trexio_file,int32_t* const n_up);

//Function reading number of molecular orbitals
trexio_exit_code trexio_read_mo_num(trexio_t* const trexio_file,int32_t* const mo_num);

//Function reading one electron integral
trexio_exit_code trexio_read_mo_1e_int_core_hamiltonian(trexio_t* const trexio_file,double* const data);

//Function reading number of non-zero two electron integral
trexio_exit_code trexio_read_mo_2e_int_eri_size(trexio_t* const trexio_file,int64_t* const n_integrals);

//Function reading the non-zero two electron integral
trexio_exit_code trexio_read_mo_2e_int_eri(trexio_t* const file,const int64_t offset_file,int64_t* const buffer_size,int32_t* const index,double* const value);

//Function reading the molecular energy
trexio_exit_code trexio_read_mo_energy(trexio_t* const file,double* const mo_energy);

// Function to calculate the Hartree-Fock energy
double calculate_hartree_fock_energy(double nuclear_repulsion_energy, double* one_e_integrals, int32_t* indices, double* two_e_values, int32_t n_occ, int32_t mo_num, int64_t n_integrals);

// Function to calculate MP2 energy correction
double calculate_mp2_energy(double* mo_energy, int32_t* indices, double* two_e_values, int32_t n_occ, int32_t mo_num, int64_t n_integrals);
#endif
