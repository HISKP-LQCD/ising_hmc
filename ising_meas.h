#ifndef EXTERN_LIBS
#define EXTERN_LIBS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "ising_aux.h"
#include "ising_evolve.h"
#endif

#ifndef ISING_MEAS
#define ISING_MEAS

double magnetization(double psi_bar, double sq_J, double h, double kappa);
double global_energy(double *phi, double *tanh_Jphi, double psi_bar, double sq_J, double h, double mass, double kappa, unsigned ns, unsigned nn);

double *measure(double *psi, double *p, unsigned *nnt, double sq_J, double h, double mass, double kappa, unsigned ns, unsigned nn, unsigned sweeps, unsigned skip, unsigned flip_freq, unsigned nmd, double dt, char *integrator, gsl_rng *r);
void write_out(double *psi, double *magn, unsigned ns, unsigned sweeps, unsigned skip, char *restart, char *save, char *out_name);
#endif
