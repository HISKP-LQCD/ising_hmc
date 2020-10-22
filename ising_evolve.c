#include "ising_evolve.h"

double hamilton(double *psi, double *phi, double *p, double psi_bar, double sq_J, double h, unsigned ns){
	unsigned i;
	double sum = 0;

	for(i = 0; i < ns; i++){
		sum += 0.5*(psi[i]*phi[i] + p[i]*p[i]);
		sum -= log(cosh(sq_J*phi[i]));
	}
	sum -= h*ns*psi_bar/sq_J;
	return sum;
}

double p_dot(double *phi, double *tanh_Jphi, double h, double sq_J, unsigned i){
	return -phi[i] + h/sq_J + sq_J*tanh_Jphi[i];
}

void leap_frog(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt){
	unsigned i, k;

	for(i = 0; i < ns; i++) psi[i] += dt/2*p[i];
	for(k = 0; k < nmd; k++){
		k_times_psi(psi, phi, mass, nnt, nnc, ns, nn);
		apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
		k_times_psi(tanh_Jphi, tanh_Jphi+ns, mass, nnt, nnc, ns, nn);
		for(i = 0; i < ns; i++){
			p[i] += dt*p_dot(phi, tanh_Jphi+ns, h, sq_J, i);
			psi[i] += dt*p[i];
		}
	}
	for(i = 0; i < ns; i++){
		psi[i] -= dt/2*p[i];
	}
	k_times_psi(psi, phi, mass, nnt, nnc, ns, nn);
	apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
}

void integrate(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, char *integrator){
	if(!strcmp(integrator, "leapfrog")){
		leap_frog(psi, p, phi, tanh_Jphi, nnt, nnc, ns, nn, h, sq_J, mass, nmd, dt);
	}else{
		printf("Integrator \"%s\" is not defined!\n", integrator);
		exit(0);
	}
}

short trajectory(double *psi, double *p, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, double *psi_bar, char *integrator, gsl_rng *r){
	unsigned i;
	double psi_bar_old = *psi_bar;
	double energy_old, energy;
	double *phi = p+ns;
	double *tanh_Jphi = p+2*ns;
	double *psi_old = p+4*ns;
	double *phi_old = p+5*ns;
	double *tanh_Jphi_old = p+6*ns;

	k_times_psi(psi, phi, mass, nnt, nnc, ns, nn);
	apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
	for(i = 0; i < ns; i++){
		p[i] = gsl_ran_ugaussian(r);
		psi_old[i] = psi[i];
		phi_old[i] = phi[i];
		tanh_Jphi_old[i] = tanh_Jphi[i];
	}

	energy_old = hamilton(psi, phi, p, *psi_bar, sq_J, h, ns);
	integrate(psi, p, phi, tanh_Jphi, nnt, nnc, ns, nn, h, sq_J, mass, nmd, dt, integrator);
	*psi_bar = average(psi, nnc, ns, nn);
	energy = hamilton(psi, phi, p, *psi_bar, sq_J, h, ns);

	if(exp(energy_old-energy) < gsl_rng_uniform(r)){
		// reject (nothing to be done if accepted)
		for(i = 0; i < ns; i++){
			psi[i] = psi_old[i];
			phi[i] = phi_old[i];
			tanh_Jphi[i] = tanh_Jphi_old[i];
		}
		*psi_bar = psi_bar_old;
		return 0;
	}else{
		return 1;
	}
}

void global_flip(double *psi, double *psi_bar, double h, double sq_J, unsigned ns, gsl_rng *r){
	const double deltaE = 2*ns*h/sq_J*(*psi_bar);
	if(exp(-deltaE) > gsl_rng_uniform(r)){
		for(unsigned i = 0; i < ns; i++) psi[i] *= -1;
		*psi_bar *= -1;
	}
}
