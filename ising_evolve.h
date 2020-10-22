#ifndef EXTERN_LIBS
#define EXTERN_LIBS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "ising_aux.h"
#endif

#ifndef ISING_ALL_GEOM
#define ISING_ALL_GEOM

double hamilton(double *psi, double *phi, double *p, double psi_bar, double sq_J, double h, unsigned ns);
double p_dot(double *phi, double *tanh_Jphi, double h, double sq_J, unsigned i);

void leap_frog(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt);

void integrate(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, char *integrator);
short trajectory(double *psi, double *p, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, double *psi_bar, char *integrator, gsl_rng *r);
void global_flip(double *psi, double *psi_bar, double h, double sq_J, unsigned ns, gsl_rng *r);
#endif
