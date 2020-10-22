#ifndef EXTERN_LIBS
#define EXTERN_LIBS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

#ifndef ISING_AUX
#define ISING_AUX

double sum_next_neighbours(double *psi, unsigned *nnt, unsigned nn, unsigned i);
double sum_next_neighbours_weighted(double *psi, unsigned *nnt, double *nnc, unsigned nn, unsigned i);
void k_times_psi(double *psi, double *phi, double mass, unsigned *nnt, double *nnc, unsigned ns, unsigned nn);
void apply_tanh_sq_J(double *phi, double *tanh_Jphi, double sq_J, unsigned ns);
double sum_vector(double* x, unsigned n);
double scalar_product(double* x, double *y, unsigned n);
double average(double *psi, double *nnc, unsigned ns, unsigned nn);
void update_fields(double *psi, double *phi, double *tanh_Jphi, double mass, double sq_J, unsigned *nnt, double *nnc, unsigned ns, unsigned nn);
#endif
