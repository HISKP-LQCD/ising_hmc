#include "ising_init.h"

gsl_rng *set_zufall(int sec){
	gsl_rng *r;
	const gsl_rng_type *T;

	gsl_rng_env_setup();
	/*langsamer, aber besser:	random8_bsd (178 MHz), random256_bsd (117 MHz),
								gfsr4 (107 MHz), mt19937 (84 MHz),
								ranlxs0 (22 MHz), ranlxs1 (15 MHz),
								ranlxd1 (8,3 MHz), ranlxd2 (4,6 MHz)*/
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, time(NULL)+sec);

	return r;
}

void print_vec(double *vec, unsigned n){
	unsigned i;
	for(i = 0; i < n; i++) printf("%g\n", vec[i]);
}

FILE *read_constants(char *file_name, double *J, double *h, double *beta, double *mass, unsigned *therm, unsigned *sweeps, double *length, unsigned *nmd, unsigned *flip_freq, int *skip, char *start, char *restart, char *save, char *out_name, char *integrator, char *geometry){
	FILE *input = fopen(file_name, "r");
	int check_sum = 0;
	if(!input){
		printf("File \"%s\" that should contain the input data, could not be opened!\n", file_name);
		exit(0);
	}
	check_sum += fscanf(input, "J=%lg\n", J);
	check_sum += fscanf(input, "h=%lg\n", h);
	check_sum += fscanf(input, "beta=%lg\n", beta);
	check_sum += fscanf(input, "mass=%lg\n", mass);
	check_sum += fscanf(input, "nr. thermalisation steps=%u\n", therm);
	check_sum += fscanf(input, "nr. trajectories=%u\n", sweeps);
	check_sum += fscanf(input, "trajectory length=%lg\n", length);
	check_sum += fscanf(input, "nr. MD-steps=%u\n", nmd);
	check_sum += fscanf(input, "inversion frequency=%u\n", flip_freq);
	check_sum += fscanf(input, "measurement frequency=%d\n", skip);
	check_sum += fscanf(input, "start type=%s\n", start);
	check_sum += fscanf(input, "start from file=%s\n", restart);
	check_sum += fscanf(input, "write config to file=%s\n", save);
	check_sum += fscanf(input, "output file name=%s\n", out_name);
	check_sum += fscanf(input, "integrator=%[^\n]s\n", integrator);
	check_sum += fscanf(input, " geometry=%s\n", geometry);
	if(check_sum != 16){
		printf("Could not read input from file \"%s\" correctly!\n", file_name);
		exit(0);
	}
	return input;
}

void set_constants(double *nnc, double *J, double *h, double beta, double *sq_J, double length, unsigned nmd, double *dt, unsigned ns, unsigned nn, double mass, double *kappa){
	*J = beta*(*J);
	*h = beta*(*h);
	*sq_J = sqrt(*J);
	*dt = length/nmd;
	if(nnc)
		*kappa = scalar_product(nnc+ns*nn, nnc+ns*(nn+1), ns)/ns;
	else
		*kappa = 1./(2*nn+mass);
}

double *initialize(double *psi, unsigned *nnt, double *nnc, double *J, double *h, double beta, double *sq_J, double length, unsigned nmd, double *dt, double mass, double *kappa, unsigned ns, unsigned nn, char *start, char *restart, unsigned therm, unsigned flip_freq, char *integrator, gsl_rng *r){
	double *container = malloc(7*ns*sizeof(double)); // {p, phi, tanh_Jphi, K*tanh_Jphi, psi_old, phi_old, tanh_Jphi_old}
	set_constants(nnc, J, h, beta, sq_J, length, nmd, dt, ns, nn, mass, kappa);
	if(strcmp(restart, "0")){
		FILE *config = fopen(restart, "r");
		if(!config){
			printf("File \"%s\" that should contain the starting config, could not be opened!\n", restart);
			exit(0);
		}
		fread(psi, sizeof(double), ns, config);
		fclose(config);
	}else{
		unsigned k;
		double psi_bar;

		if(!strcmp(start, "hot")){
			for(k = 0; k < ns; k++) psi[k] = gsl_ran_gaussian_ziggurat(r, *sq_J);
		}else if(!strcmp(start, "cold")){
			double psi_0 = (*h >= 0? 1: -1) * (*sq_J) + (*h)*(*kappa)/(*sq_J);
			for(k = 0; k < ns; k++) psi[k] = psi_0;
		}else if(!strcmp(start, "zero")){
			for(k = 0; k < ns; k++) psi[k] = 0;
		}else{
			printf("Start type \"%s\" is not defined!\n", start);
			exit(0);
		}

		psi_bar = average(psi, nnc, ns, nn);

		for(k = 0; k < therm; k++){
			trajectory(psi, container, nnt, nnc, ns, nn, *h, *sq_J, mass, nmd, *dt, &psi_bar, integrator, r);
			if(k%flip_freq == 0) global_flip(psi, &psi_bar, *h, *sq_J, ns, r);
		}
	}
	return container;
}

void free_all(double *psi, double *p, unsigned *nnt, double *nnc, double *magn){
	free(psi);
	free(p);
	free(nnt);
	free(nnc);
	free(magn);
}
