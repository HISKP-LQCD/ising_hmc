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
#endif

#ifndef ISING_ALL_GEOM
#define ISING_ALL_GEOM

void sample_fourier_phi(double *phi, double complex *xc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, gsl_rng *r, const fftw_plan *fft);
void sample_fourier_momenta(double *p, double complex *pc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, gsl_rng *r, const fftw_plan *fft);
double energy_fourier_phi(double complex *xc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, const fftw_plan *fft);
double energy_fourier_momenta(double complex *pc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, const fftw_plan *fft);
void evolve_bosons(double complex *xc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, const fftw_plan *fft, double dt);
double hamilton(double complex *psi, double *phi, double *p, double mass, double sq_J, double h, unsigned ns, unsigned nn, unsigned nl, unsigned dim, const fftw_plan *fft);
double p_dot(double *tanh_Jphi, double h, double sq_J, unsigned i);
void leap_frog(double complex *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, unsigned nl, unsigned dim, double h, double sq_J, double mass, unsigned nmd, double dt, const fftw_plan *fft);
void integrate(double complex *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, unsigned nl, unsigned dim, double h, double sq_J, double mass, unsigned nmd, double dt, char *integrator, const fftw_plan *fft);
short trajectory(double complex *psi, double *p, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, unsigned nl, unsigned dim, double h, double sq_J, double mass, unsigned nmd, double dt, double *psi_bar, char *integrator, gsl_rng *r, const fftw_plan *fft);
void global_flip(double *psi, double *psi_bar, double h, double sq_J, unsigned ns, gsl_rng *r);
#endif
