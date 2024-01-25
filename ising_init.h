#ifndef EXTERN_LIBS
#define EXTERN_LIBS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <fftw3.h>

#include "ising_lattice.h"
#include "ising_aux.h"
#include "ising_evolve.h"
#endif

#ifndef ISING_INIT
#define ISING_INIT

gsl_rng *set_zufall(int sec);

FILE *read_constants(char *file_name, double *J, double *h, double *beta, double *mass, unsigned *therm, unsigned *sweeps, double *length, unsigned *nmd, unsigned *flip_freq, int *skip, char *start, char *restart, char *save, char *out_name, char *integrator, char *geometry);
void set_constants(double *nnc, double *J, double *h, double beta, double *sq_J, double length, unsigned nmd, double *dt, unsigned ns, unsigned nn, double mass, double *kappa);

double *initialize(double complex *psi, unsigned *nnt, double *nnc, double *J, double *h, double beta, double *sq_J, double length, unsigned nmd, double *dt, double mass, double *kappa, unsigned ns, unsigned nn, unsigned nl, unsigned dim, char *start, char *restart, unsigned therm, unsigned flip_freq, char *integrator, gsl_rng *r, fftw_plan *fft);
void free_all(double complex *psi, double *p, unsigned *nnt, double *nnc, double *magn, fftw_plan *fft);
#endif
