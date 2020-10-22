#include "ising_meas.h"

double magnetization(double psi_bar, double sq_J, double h, double kappa){
	return psi_bar/sq_J - h*kappa/(sq_J*sq_J);
}

double global_energy(double *phi, double *tanh_Jphi, double psi_bar, double sq_J, double h, double mass, double kappa, unsigned ns, unsigned nn){
	double sum = sq_J*sq_J*(nn+mass)/2
				+ h*h*kappa/(2*sq_J*sq_J)
				- h*psi_bar/(2*sq_J)
				- sq_J*scalar_product(phi, tanh_Jphi, ns)/(2*ns);
	return sum;
}

double *measure(double *psi, double *p, unsigned *nnt, double *nnc, double sq_J, double h, double mass, double kappa, unsigned ns, unsigned nn, unsigned sweeps, unsigned skip, unsigned flip_freq, unsigned nmd, double dt, char *integrator, gsl_rng *r){
	unsigned nr_meas = sweeps/skip;
	unsigned i, k;
	short acc;
	double psi_bar = average(psi, nnc, ns, nn);
	double *magn = malloc(3*nr_meas * sizeof(double));
	double *energy = magn+nr_meas;
	double *accepted = magn+2*nr_meas;

	for(i = 0; i < sweeps; i++){
		acc = trajectory(psi, p, nnt, nnc, ns, nn, h, sq_J, mass, nmd, dt, &psi_bar, integrator, r);
		if(i%flip_freq == 0) global_flip(psi, &psi_bar, h, sq_J, ns, r);
		if(i%skip == 0){
			k = i/skip;
			magn[k] = magnetization(psi_bar, sq_J, h, kappa);
			energy[k] = global_energy(p+ns, p+2*ns, psi_bar, sq_J, h, mass, kappa, ns, nn);
			accepted[k] = acc;
		}
	}

	return magn;
}

void write_out(double *psi, double *magn, unsigned ns, unsigned sweeps, unsigned skip, char *restart, char *save, char *out_name){
	const unsigned nr_meas = sweeps/skip;
	unsigned i;
	double *energy = magn+nr_meas;
	double *accepted = magn+2*nr_meas;
	FILE *out;

	if(strcmp(save, "0")){
		FILE *config = fopen(save, "w");
		fwrite(psi, sizeof(double), ns, config);
		fclose(config);
	}

	if(strcmp(restart, "0")){
		out = fopen(out_name, "a");
	}else{
		out = fopen(out_name, "w");
	}
	for(i = 0; i < nr_meas; i++){
		const double m = magn[i];
		fprintf(out, "%.15g\t%.15g\t%.15g\t%.15g\t%.0f\n", m, fabs(m), m*m, energy[i], accepted[i]);
	}
	fclose(out);
}
