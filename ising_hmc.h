#ifndef EXTERN_LIBS
#define EXTERN_LIBS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <fftw3.h>

#include "ising_init.h"
#include "ising_lattice.h"
#include "ising_aux.h"
#include "ising_evolve.h"
#include "ising_meas.h"
#endif

#ifndef ISING_HMC
#define ISING_HMC

void display_help();
#endif
