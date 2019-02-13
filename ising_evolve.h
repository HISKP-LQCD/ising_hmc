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
double p_dot(double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double weight, unsigned i);
double p_dot_orthogonal(double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double sq_J, double weight, unsigned i, double shift);
double p_dot_parallel(double psi0, double delta_psi, double weight, double h, double sq_J, double *tanh_averages, unsigned truncation_order);

void evolve_Omelyan_4th_parallel(double *psi, double *p, double weight, double h, double sq_J, double *tanh_averages, unsigned nmd, double dt, unsigned truncation_order);
void evolve_Omelyan_4th_split(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, unsigned sub_nmd, double sub_dt, unsigned truncation_order);
void evolve_Omelyan_4th(double *psi, double *pi, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt);

void leap_frog_parallel(double *psi, double *d_psi0, double *p, double weight, double h, double sq_J, double *tanh_averages, unsigned nmd, double dt, unsigned truncation_order, int up);
void leap_frog_split(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, unsigned sub_nmd, double sub_dt, unsigned truncation_order);
void leap_frog(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt);

void integrate(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, char *integrator);
short trajectory(double *psi, double *p, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, double *psi_bar, char *integrator, gsl_rng *r);
void global_flip(double *psi, double *psi_bar, double h, double sq_J, unsigned ns, gsl_rng *r);
#endif
