#include "ising_aux.h"

double sum_next_neighbours(double *psi, unsigned *nnt, unsigned nn, unsigned i){
	unsigned k;
	double sum_nn = 0;

	nnt += i*nn;
	for(k = 0; k < nn; k++) sum_nn += psi[nnt[k]];
	return sum_nn;
}

double sum_next_neighbours_weighted(double *psi, unsigned *nnt, double *nnc, unsigned nn, unsigned i){
	unsigned k;
	double sum_nn = 0;

	nnt += i*nn;
	nnc += i*nn;
	for(k = 0; k < nn; k++) sum_nn += nnc[k]*psi[nnt[k]];
	return sum_nn;
}

void k_times_psi(double *psi, double *phi, double mass, unsigned *nnt, double *nnc, unsigned ns, unsigned nn){
	unsigned i;
	const double weight = nn+mass;

	if(nnc)
		for(i = 0; i < ns; i++) phi[i] = weight*psi[i] + sum_next_neighbours_weighted(psi, nnt, nnc, nn, i);
	else
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

double average(double *psi, double *nnc, unsigned ns, unsigned nn){
	if(nnc)
		return scalar_product(psi, nnc+ns*nn, ns)/ns;
	else
		return sum_vector(psi, ns)/ns;
}

void update_fields(double *psi, double *phi, double *tanh_Jphi, double mass, double sq_J, unsigned *nnt, double *nnc, unsigned ns, unsigned nn){
	k_times_psi(psi, phi, mass, nnt, nnc, ns, nn);
	apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
}
