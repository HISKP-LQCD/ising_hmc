#include "ising_aux.h"

double sum_next_neighbours(double *psi, unsigned *nnt, unsigned nn, unsigned i){
	unsigned k;
	double sum_nn = 0;

	nnt += i*nn;
	for(k = 0; k < nn; k++) sum_nn += psi[nnt[k]];
	return sum_nn;
}

void k_times_psi(double *psi, double *phi, double mass, unsigned *nnt, unsigned ns, unsigned nn){
	unsigned i;
	double weight = nn+mass;

	for(i = 0; i < ns; i++) phi[i] = weight*psi[i] + sum_next_neighbours(psi, nnt, nn, i);
}

void apply_tanh_sq_J(double *phi, double *tanh_Jphi, double sq_J, unsigned ns){
	unsigned i;
	for(i = 0; i < ns; i++) tanh_Jphi[i] = tanh(sq_J * phi[i]);
}

double sum_vector(double* x, unsigned n){
	unsigned i;
	double sum = 0;
	for(i = 0; i < n; i++) sum += x[i];
	return sum;
}

double scalar_product(double* x, double *y, unsigned n){
	unsigned i;
	double sum = 0;
	for(i = 0; i < n; i++) sum += x[i]*y[i];
	return sum;
}

void expand_tanh(double *tanh_Jphi, double *tanh_averages, unsigned ns, unsigned truncation_order){
	double th_phi_toi;
	unsigned i, k;
	for(k = 0; k < truncation_order; k++) tanh_averages[k] = 0;
	for(i = 0; i < ns; i++){
		const double th_phi = tanh_Jphi[i];
		const double th_sq_minus1 = th_phi*th_phi-1;
		th_phi_toi = 1;
		tanh_averages[0] += th_phi;
		tanh_averages[1] += th_sq_minus1;
		for(k = 2; k < truncation_order; k++){
			th_phi_toi *= th_phi;
			tanh_averages[k] += th_phi_toi*th_sq_minus1;
		}
	}
	for(k = 0; k < truncation_order; k++) tanh_averages[k] /= ns;
}

void update_fields(double *psi, double *phi, double *tanh_Jphi, double mass, double sq_J, unsigned *nnt, unsigned ns, unsigned nn){
	k_times_psi(psi, phi, mass, nnt, ns, nn);
	apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
}

void update_fields_split(double *psi, double *phi, double *tanh_Jphi, double *tanh_averages, double mass, double sq_J, unsigned *nnt, unsigned ns, unsigned nn, unsigned truncation_order){
	k_times_psi(psi, phi, mass, nnt, ns, nn);
	apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
	expand_tanh(tanh_Jphi, tanh_averages, ns, truncation_order);
}
