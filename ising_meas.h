#ifndef EXTERN_LIBS
#define EXTERN_LIBS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>

#include "ising_aux.h"
#include "ising_evolve.h"
#endif

#ifndef ISING_MEAS
#define ISING_MEAS

FILE *file_state(char *out_name, char *restart, int skip);
void print_state(FILE *out, double *psi, double *nnc, double sq_J, double h, double kappa, unsigned ns, unsigned nn, int skip);
double magnetization(double psi_bar, double sq_J, double h, double kappa);
double global_energy(double *phi, double *tanh_Jphi, double psi_bar, double sq_J, double h, double mass, double kappa, unsigned ns, unsigned nn);

double *measure(double complex *psi, double *p, unsigned *nnt, double *nnc, double sq_J, double h, double mass, double kappa, unsigned ns, unsigned nn, unsigned nl, unsigned dim, unsigned sweeps, int skip, unsigned flip_freq, unsigned nmd, double dt, char *integrator, char *restart, char *out_name, gsl_rng *r, const fftw_plan *fft);
void write_out(double *psi, double *magn, unsigned ns, unsigned sweeps, int skip, char *restart, char *save, char *out_name);
#endif
