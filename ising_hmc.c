// make
// ./ising_hmc <input-file>
//
// Author: Johann Ostmeyer

#include "ising_hmc.h"

void display_help(){
	printf("You have to supply at least one file with input data!\n");
	printf("The input file has to have the following format:\n\n");
	printf("J=<coupling>\nh=<external field>\nbeta=<inverse temperature>\nnr. thermalisation steps=<updates performed before measurement>\nnr. sweeps=<length of measurement>\nmeasurement frequency=<skip sweeps for less autocorrelation, negative for saving the state as well>\nstart type=<hot/cold/zero>\nstart from file=<file with thermalized config, or 0 for fresh start>\nwrite config to file=<file to store the last config, or 0 if not required>\noutput file name=<file for measurements, they will be saved in columns: m | abs(m) | m^2 | beta*E>\ngeometry=<cubic/rectangular/triangular/alltoall/generic>\n");
	printf("<Further geometry-specific information>\n\n");
	printf("See function construct_lattice() and called constructors in code for more details.\n");
	printf("Exiting.\n");
}

int main(int argc, char **argv){
	gsl_rng *r=set_zufall(0);
	double J, h, beta;
	double mass, length, sq_J, dt, kappa;
	unsigned therm, sweeps;
	int skip;
	unsigned nmd, flip_freq;
	unsigned ns, nn;
	unsigned i;
	char start[100], restart[500], save[500], out_name[500], integrator[500], geometry[200];
	double *psi, *p;
	unsigned *nnt;
	double *nnc = NULL;
	FILE *input;
	double *magn;

	if(argc < 2){
		display_help();
		return 0;
	}

	for(i = 1; i < argc; i++){
		input = read_constants(argv[i], &J, &h, &beta, &mass, &therm, &sweeps, &length, &nmd, &flip_freq, &skip, start, restart, save, out_name, integrator, geometry);
		printf("Read constants.\n");
		nnt = construct_lattice(geometry, input, &ns, &nn, &psi, &nnc);
		printf("Constructed lattice.\n");
		p = initialize(psi, nnt, nnc, &J, &h, beta, &sq_J, length, nmd, &dt, mass, &kappa, ns, nn, start, restart, therm, flip_freq, integrator, r);
		printf("Initialized and thermalized.\n");
		magn = measure(psi, p, nnt, nnc, sq_J, h, mass, kappa, ns, nn, sweeps, skip, flip_freq, nmd, dt, integrator, restart, out_name, r);
		printf("Measured.\n");
		write_out(psi, magn, ns, sweeps, skip, restart, save, out_name);
		printf("Wrote out results.\n");
		free_all(psi, p, nnt, nnc, magn);
		printf("Finished %d calculations.\n\n", i);
	}

	gsl_rng_free(r);

	return 0;
}
