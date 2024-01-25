#include "ising_evolve.h"

void sample_fourier_phi(double *phi, double complex *xc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, gsl_rng *r, const fftw_plan *fft){
	const unsigned loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	const double matsubara = 2 * M_PI / nl;

	for(unsigned i = 0; i < ns; i++)
		phi[i] = gsl_ran_ugaussian(r);

	fftw_execute(fft[0]);

	for(unsigned i = 0; i < compl_dim; i++){
		double freq = nn + mass;
		unsigned ind = i / loc_dim, pos = i % loc_dim;

		for(unsigned k = 0; k < dim; k++){
			freq += 2 * cos(matsubara * pos);
			pos = ind % nl;
			ind /= nl;
		}

		const double ampl = sqrt(freq) / ns; // factor ns because fft accumulates it from there-back trafo

		xc[i] *= ampl;
	}

	fftw_execute(fft[2]);
}

void sample_fourier_momenta(double *p, double complex *pc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, gsl_rng *r, const fftw_plan *fft){
	const unsigned loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	const double matsubara = 2 * M_PI / nl;

	for(unsigned i = 0; i < ns; i++)
		p[i] = gsl_ran_ugaussian(r);

	fftw_execute(fft[1]);

	for(unsigned i = 0; i < compl_dim; i++){
		double freq = nn + mass;
		unsigned ind = i / loc_dim, pos = i % loc_dim;

		for(unsigned k = 0; k < dim; k++){
			freq += 2 * cos(matsubara * pos);
			pos = ind % nl;
			ind /= nl;
		}

		const double ampl = 1. / sqrt(freq) / ns; // factor ns because fft accumulates it from there-back trafo

		pc[i] *= ampl;
	}

	fftw_execute(fft[3]);
}

double energy_fourier_phi(double complex *xc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, const fftw_plan *fft){
	const unsigned loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	const double matsubara = 2 * M_PI / nl;
	double en = 0;

	fftw_execute(fft[0]);

	for(unsigned i = 0; i < compl_dim; i++){
		unsigned ind = i / loc_dim, pos = i % loc_dim;
		double freq = nn + mass, weight = .5;

		// elements t=1...Nt/2-1 occur twice (as complex conj pairs), but are stored only once
		if(pos > 0 && pos < (nl+1)/2) weight *= 2;

		for(unsigned k = 0; k < dim; k++){
			freq += 2 * cos(matsubara * pos);
			pos = ind % nl;
			ind /= nl;
		}

		weight *= 1. / freq / ns; // factor ns because fft accumulates it from there-back trafo

		en += norm(xc[i]) * weight;
	}

	return en;
}

double energy_fourier_momenta(double complex *pc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, const fftw_plan *fft){
	const unsigned loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	const double matsubara = 2 * M_PI / nl;
	double en = 0;

	fftw_execute(fft[1]);

	for(unsigned i = 0; i < compl_dim; i++){
		unsigned ind = i / loc_dim, pos = i % loc_dim;
		double freq = nn + mass, weight = .5;

		// elements t=1...Nt/2-1 occur twice (as complex conj pairs), but are stored only once
		if(pos > 0 && pos < (nl+1)/2) weight *= 2;

		for(unsigned k = 0; k < dim; k++){
			freq += 2 * cos(matsubara * pos);
			pos = ind % nl;
			ind /= nl;
		}

		weight *= freq / ns; // factor ns because fft accumulates it from there-back trafo

		en += norm(pc[i]) * weight;
	}

	return en;
}

void evolve_bosons(double complex *xc, double mass, unsigned ns, unsigned nn, unsigned nl, unsigned dim, const fftw_plan *fft, double dt){
	const unsigned loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	const double matsubara = 2 * M_PI / nl;
	const double c = cos(dt) / ns, s = sin(dt) / ns; // factor ns because fft accumulates it from there-back trafo
	double complex *pc = xc + compl_dim;

	fftw_execute(fft[0]);
	fftw_execute(fft[1]);

	for(unsigned i = 0; i < compl_dim; i++){
		double freq = nn + mass;
		unsigned ind = i / loc_dim, pos = i % loc_dim;

		for(unsigned k = 0; k < dim; k++){
			freq += 2 * cos(matsubara * pos);
			pos = ind % nl;
			ind /= nl;
		}

		const double sw = s * freq, s1w = s / freq;
		const double complex x0 = xc[i], p0 = pc[i];
		xc[i] = c * x0 + sw  * p0;
		pc[i] = c * p0 - s1w * x0;
	}

	fftw_execute(fft[2]);
	fftw_execute(fft[3]);
}

double hamilton(double complex *psi, double *phi, double *p, double mass, double sq_J, double h, unsigned ns, unsigned nn, unsigned nl, unsigned dim, const fftw_plan *fft){
	const unsigned loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	unsigned i;
	double sum = 0;

	sum += energy_fourier_momenta(psi + compl_dim, mass, ns, nn, nl, dim, fft);
	//printf("E_p = %g ,\t", sum);
	sum += energy_fourier_phi(psi, mass, ns, nn, nl, dim, fft);
	//printf("+E_phi = %g ,\t", sum);

	for(i = 0; i < ns; i++){
		sum -= log(cosh(sq_J*phi[i] + h));
	}
	//printf("+E_m = %g\n", sum);
	return sum;
}

double p_dot(double *tanh_Jphi, double h, double sq_J, unsigned i){
	return sq_J*tanh_Jphi[i];
}

void leap_frog(double complex *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, unsigned nl, unsigned dim, double h, double sq_J, double mass, unsigned nmd, double dt, const fftw_plan *fft){
	unsigned i, k;

	if(!nmd) dt *= 2;
	evolve_bosons(psi, mass, ns, nn, nl, dim, fft, .5*dt);

	for(k = 1; k <= nmd; k++){
		apply_tanh_sq_J(phi, tanh_Jphi, sq_J, h, ns);
		for(i = 0; i < ns; i++){
			p[i] += dt*p_dot(tanh_Jphi, h, sq_J, i);
		}

		if(k == nmd) dt *= .5; // only half step in the end
		evolve_bosons(psi, mass, ns, nn, nl, dim, fft, dt);
	}

	apply_tanh_sq_J(phi, tanh_Jphi, sq_J, h, ns);
}

void integrate(double complex *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, unsigned nl, unsigned dim, double h, double sq_J, double mass, unsigned nmd, double dt, char *integrator, const fftw_plan *fft){
	if(!strcmp(integrator, "leapfrog")){
		leap_frog(psi, p, phi, tanh_Jphi, nnt, nnc, ns, nn, nl, dim, h, sq_J, mass, nmd, dt, fft);
	}else{
		printf("Integrator \"%s\" is not defined!\n", integrator);
		exit(0);
	}
}

short trajectory(double complex *psi, double *p, unsigned *nnt, double *nnc, unsigned ns, unsigned nn, unsigned nl, unsigned dim, double h, double sq_J, double mass, unsigned nmd, double dt, double *psi_bar, char *integrator, gsl_rng *r, const fftw_plan *fft){
	const unsigned loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	unsigned i;
	double psi_bar_old = psi_bar? *psi_bar: 0;
	double energy_old, energy;
	double *phi = p+ns;
	double *tanh_Jphi = p+2*ns;
	double *phi_old = p+5*ns;
	double *tanh_Jphi_old = p+6*ns;

	apply_tanh_sq_J(phi, tanh_Jphi, sq_J, h, ns);
	for(i = 0; i < ns; i++){
		phi_old[i] = phi[i];
		tanh_Jphi_old[i] = tanh_Jphi[i];
	}
	sample_fourier_momenta(p, psi + compl_dim, mass, ns, nn, nl, dim, r, fft);

	energy_old = hamilton(psi, phi, p, mass, sq_J, h, ns, nn, nl, dim, fft);
	integrate(psi, p, phi, tanh_Jphi, nnt, nnc, ns, nn, nl, dim, h, sq_J, mass, nmd, dt, integrator, fft);
	if(psi_bar) *psi_bar = average(tanh_Jphi, nnc, ns, nn);
	energy = hamilton(psi, phi, p, mass, sq_J, h, ns, nn, nl, dim, fft);

	printf("%g\t%g\t%g\t%g\n", energy_old, energy, energy-energy_old, exp(energy_old-energy));

	if(exp(energy_old-energy) < gsl_rng_uniform(r)){
		// reject (nothing to be done if accepted)
		for(i = 0; i < ns; i++){
			phi[i] = phi_old[i];
			tanh_Jphi[i] = tanh_Jphi_old[i];
		}
		if(psi_bar) *psi_bar = psi_bar_old;
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
