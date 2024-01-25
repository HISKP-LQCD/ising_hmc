#ifndef EXTERN_LIBS
#define EXTERN_LIBS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#endif

#ifndef ISING_LATTICE
#define ISING_LATTICE

unsigned *construct_lattice(char *geometry, FILE *input, unsigned *ns, unsigned *nn, unsigned *nl, unsigned *dim, double **s, double **nnc);
unsigned *construct_cubic(FILE *input, unsigned *ns, unsigned *nn, unsigned *l, unsigned *dim, double **s);
unsigned *construct_rectangular(FILE *input, unsigned *ns, unsigned *nn, double **s);
unsigned *construct_triangular(FILE *input, unsigned *ns, unsigned *nn, double **s);
unsigned *construct_alltoall(FILE *input, unsigned *ns, unsigned *nn, double **s);
unsigned *construct_generic(FILE *input, unsigned *ns, unsigned *nn, double **s, double **nnc);
#endif
