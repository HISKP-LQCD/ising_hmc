#include "ising_meas.h"

FILE *file_state(char *out_name, char *restart, int skip){
	char out_state[500];
	FILE *out = NULL;

	if(skip < 0){
		sprintf(out_state, "%s.state", out_name);
		if(strcmp(restart, "0")){
			out = fopen(out_state, "a");
		}else{
			out = fopen(out_state, "w");
		}
	}
	return out;
}

void print_state(FILE *out, double *psi, double *nnc, double sq_J, double h, double kappa, unsigned ns, unsigned nn, int skip){
	if(skip < 0){
		const double inv_sq_J = 1./sq_J;
		if(nnc){
			const double shift = h*inv_sq_J*inv_sq_J;
			nnc += ns*(nn+1);
			for(unsigned i = 0; i < ns; i++) fprintf(out, "%.15g\t", psi[i]*inv_sq_J - nnc[i]*shift);
		}else{
			const double shift = h*kappa/(sq_J*sq_J);
			for(unsigned i = 0; i < ns; i++) fprintf(out, "%.15g\t", psi[i]*inv_sq_J - shift);
		}
		fprintf(out, "\n");
	}
}

double magnetization(double psi_bar, double sq_J, double h, double kappa){
	// deprecated
	return psi_bar/sq_J - h*kappa/(sq_J*sq_J);
}

double global_energy(double *phi, double *tanh_Jphi, double psi_bar, double sq_J, double h, double mass, double kappa, unsigned ns, unsigned nn){
	double sum = sq_J*sq_J*(nn+mass)/2
				- sq_J*scalar_product(phi, tanh_Jphi, ns)/(2*ns)
	            - h*psi_bar/2;
	return sum;
}

double *measure(double complex *psi, double *p, unsigned *nnt, double *nnc, double sq_J, double h, double mass, double kappa, unsigned ns, unsigned nn, unsigned nl, unsigned dim, unsigned sweeps, int skip, unsigned flip_freq, unsigned nmd, double dt, char *integrator, char *restart, char *out_name, gsl_rng *r, const fftw_plan *fft){
	const unsigned skip_freq = abs(skip);
	const unsigned nr_meas = sweeps/skip_freq;
	unsigned i, k;
	short acc;
	double psi_bar = average(p+2*ns, nnc, ns, nn);
	double *magn = malloc(3*nr_meas * sizeof(double));
	double *energy = magn+nr_meas;
	double *accepted = magn+2*nr_meas;
	FILE *out_state = file_state(out_name, restart, skip);

	for(i = 0; i < sweeps; i++){
		acc = trajectory(psi, p, nnt, nnc, ns, nn, nl, dim, h, sq_J, mass, nmd, dt, &psi_bar, integrator, r, fft);
		//if(i%flip_freq == 0) global_flip(psi, &psi_bar, h, sq_J, ns, r);
		if(i%skip_freq == 0){
			k = i/skip_freq;
			magn[k] = psi_bar;
			energy[k] = global_energy(p+ns, p+2*ns, psi_bar, sq_J, h, mass, kappa, ns, nn);
			accepted[k] = acc;

			print_state(out_state, p+ns, nnc, sq_J, h, kappa, ns, nn, skip);
		}
	}

	if(out_state) fclose(out_state);

	return magn;
}

void write_out(double *phi, double *magn, unsigned ns, unsigned sweeps, int skip, char *restart, char *save, char *out_name){
	const unsigned skip_freq = abs(skip);
	const unsigned nr_meas = sweeps/skip_freq;
	unsigned i;
	double *energy = magn+nr_meas;
	double *accepted = magn+2*nr_meas;
	FILE *out;

	if(strcmp(save, "0")){
		FILE *config = fopen(save, "w");
		fwrite(phi, sizeof(double), ns, config);
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
