#include "ising_evolve.h"

double hamilton(double *psi, double *phi, double *p, double psi_bar, double sq_J, double h, unsigned ns){
	unsigned i;
	double sum = 0;

	for(i = 0; i < ns; i++){
		sum += 0.5*(psi[i]*phi[i] + p[i]*p[i]);
		sum -= log(cosh(sq_J*phi[i]));
	}
	return sum - ns*h*psi_bar/sq_J;
}

double p_dot(double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double weight, unsigned i){
	return -phi[i] + h/sq_J + sq_J*(weight*tanh_Jphi[i] + sum_next_neighbours(tanh_Jphi, nnt, nn, i));
}

double p_dot_orthogonal(double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double sq_J, double weight, unsigned i, double shift){
	return -phi[i] + sq_J*(weight*tanh_Jphi[i] + sum_next_neighbours(tanh_Jphi, nnt, nn, i)) - shift;
}

double p_dot_parallel(double psi0, double delta_psi, double weight, double h, double sq_J, double *tanh_averages, unsigned truncation_order){
	const double th_psi0 = -tanh(sq_J*weight*delta_psi);
	double th_psi0_toi = 1;
	double sum_expansion = tanh_averages[0];
	unsigned i;
	for(i = 1; i < truncation_order; i++){
		th_psi0_toi *= th_psi0;
		sum_expansion += th_psi0_toi*tanh_averages[i];
	}
	return weight*(-(psi0 + delta_psi) + sq_J*sum_expansion) + h/sq_J;
}

void evolve_Omelyan_4th_parallel(double *psi, double *p, double weight, double h, double sq_J, double *tanh_averages, unsigned nmd, double dt, unsigned truncation_order){
	const double c1 = 0.08398315262876693*dt; // vartheta
	const double d1 = 0.2539785108410595*dt; // rho
	const double c2 = 0.6822365335719091*dt; // lambda
	const double d2 = -0.03230286765269967*dt; // theta
	const double c3 = 0.5*dt-(c1+c2);
	const double d3 = dt-2*(d1+d2);
	const double c12 = 2*c1;
	const double psi0 = *psi;
	double d_psi = 0, p0 = *p;
	unsigned k;

	p0 += c1 * p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order);
	for(k = 0; k < nmd; k++){
		d_psi += d1 * p0;
		p0 += c2 * p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order);
		d_psi += d2 * p0;
		p0 += c3 * p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order);
		d_psi += d3 * p0;
		p0 += c3 * p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order);
		d_psi += d2 * p0;
		p0 += c2 * p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order);
		d_psi += d1 * p0;
		p0 += c12 * p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order);
	}
	p0 -= c1 * p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order);

	*psi = psi0 + d_psi;
	*p = p0;
}

void evolve_Omelyan_4th_split(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, unsigned sub_nmd, double sub_dt, unsigned truncation_order){
	const double c1 = 0.08398315262876693*dt; // vartheta
	const double d1 = 0.2539785108410595*dt; // rho
	const double c2 = 0.6822365335719091*dt; // lambda
	const double d2 = -0.03230286765269967*dt; // theta
	const double c3 = 0.5*dt-(c1+c2);
	const double d3 = dt-2*(d1+d2);
	const double c12 = 2*c1;

	const double sub1 = 0.08398315262876693*sub_dt;
	const double sub2 = 0.1699953582122926*sub_dt;
	const double sub3 = 0.5122411753596166*sub_dt;
	const double sub4 = -0.5445440430123163*sub_dt;
	const double sub5 = 0.2783243568116402*sub_dt;

	const double weight = nn+mass, weight_parallel=2*nn+mass;
	double psi0=sum_vector(psi,ns)/ns, d_psi, p0=sum_vector(p,ns)/ns, shift;
	double *tanh_averages = malloc(truncation_order*sizeof(double));
	unsigned i, k;

	expand_tanh(tanh_Jphi, tanh_averages, ns, truncation_order);
	shift = weight_parallel*(-psi0 + sq_J*tanh_averages[0]);

	// t = 0
	for(i = 0; i < ns; i++){
		p[i] += c1 * p_dot_orthogonal(phi, tanh_Jphi, nnt, ns, nn, sq_J, weight, i, shift) - p0;
	}
	for(k = 0; k < nmd; k++){
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub1, truncation_order);
		for(i = 0; i < ns; i++) psi[i] += d1 * (p[i] + p0);
		// t = rho * dt
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub2, truncation_order);
		d_psi = psi0 - sum_vector(psi, ns)/ns;
		for(i = 0; i < ns; i++) psi[i] += d_psi;
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		shift = weight_parallel*(-psi0 + sq_J*tanh_averages[0]);
		for(i = 0; i < ns; i++) p[i] += c2 * p_dot_orthogonal(phi, tanh_Jphi, nnt, ns, nn, sq_J, weight, i, shift);
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub3, truncation_order);
		for(i = 0; i < ns; i++) psi[i] += d2 * (p[i] + p0);
		// t = (theta + rho) * dt
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub4, truncation_order);
		d_psi = psi0 - sum_vector(psi, ns)/ns;
		for(i = 0; i < ns; i++) psi[i] += d_psi;
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		shift = weight_parallel*(-psi0 + sq_J*tanh_averages[0]);
		for(i = 0; i < ns; i++) p[i] += c3 * p_dot_orthogonal(phi, tanh_Jphi, nnt, ns, nn, sq_J, weight, i, shift);
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub5, truncation_order);
		for(i = 0; i < ns; i++) psi[i] += d3 * (p[i] + p0);
		// t = (1 - 2(theta + rho)) * dt
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub5, truncation_order);
		d_psi = psi0 - sum_vector(psi, ns)/ns;
		for(i = 0; i < ns; i++) psi[i] += d_psi;
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		shift = weight_parallel*(-psi0 + sq_J*tanh_averages[0]);
		for(i = 0; i < ns; i++) p[i] += c3 * p_dot_orthogonal(phi, tanh_Jphi, nnt, ns, nn, sq_J, weight, i, shift);
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub4, truncation_order);
		for(i = 0; i < ns; i++) psi[i] += d2 * (p[i] + p0);
		// t = (1 - (theta + rho)) * dt
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub3, truncation_order);
		d_psi = psi0 - sum_vector(psi, ns)/ns;
		for(i = 0; i < ns; i++) psi[i] += d_psi;
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		shift = weight_parallel*(-psi0 + sq_J*tanh_averages[0]);
		for(i = 0; i < ns; i++) p[i] += c2 * p_dot_orthogonal(phi, tanh_Jphi, nnt, ns, nn, sq_J, weight, i, shift);
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub2, truncation_order);
		for(i = 0; i < ns; i++) psi[i] += d1 * (p[i] + p0);
		// t = (1 - rho) * dt
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		evolve_Omelyan_4th_parallel(&psi0, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd, sub1, truncation_order);
		d_psi = psi0 - sum_vector(psi, ns)/ns;
		for(i = 0; i < ns; i++) psi[i] += d_psi;
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		shift = weight_parallel*(-psi0 + sq_J*tanh_averages[0]);
		for(i = 0; i < ns; i++) p[i] += c12 * p_dot_orthogonal(phi, tanh_Jphi, nnt, ns, nn, sq_J, weight, i, shift);
	}
	for(i = 0; i < ns; i++){
		p[i] += -c1 * p_dot_orthogonal(phi, tanh_Jphi, nnt, ns, nn, sq_J, weight, i, shift) + p0;
	}

	free(tanh_averages);
}

void evolve_Omelyan_4th(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt){
	const double c1 = 0.08398315262876693*dt; // vartheta
	const double d1 = 0.2539785108410595*dt; // rho
	const double c2 = 0.6822365335719091*dt; // lambda
	const double d2 = -0.03230286765269967*dt; // theta
	const double c3 = 0.5*dt-(c1+c2);
	const double d3 = dt-2*(d1+d2);
	const double c12 = 2*c1;
	const double weight = nn+mass;
	unsigned i, k;

	for(i = 0; i < ns; i++){
		p[i] += c1 * p_dot(phi, tanh_Jphi, nnt, ns, nn, h, sq_J, weight, i);
	}
	for(k = 0; k < nmd; k++){
		for(i = 0; i < ns; i++){
			psi[i] += d1 * p[i];
		}
		k_times_psi(psi, phi, mass, nnt, ns, nn);
		apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
		for(i = 0; i < ns; i++){
			p[i] += c2 * p_dot(phi, tanh_Jphi, nnt, ns, nn, h, sq_J, weight, i);
			psi[i] += d2 * p[i];
		}
		k_times_psi(psi, phi, mass, nnt, ns, nn);
		apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
		for(i = 0; i < ns; i++){
			p[i] += c3 * p_dot(phi, tanh_Jphi, nnt, ns, nn, h, sq_J, weight, i);
			psi[i] += d3 * p[i];
		}
		k_times_psi(psi, phi, mass, nnt, ns, nn);
		apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
		for(i = 0; i < ns; i++){
			p[i] += c3 * p_dot(phi, tanh_Jphi, nnt, ns, nn, h, sq_J, weight, i);
			psi[i] += d2 * p[i];
		}
		k_times_psi(psi, phi, mass, nnt, ns, nn);
		apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
		for(i = 0; i < ns; i++){
			p[i] += c2 * p_dot(phi, tanh_Jphi, nnt, ns, nn, h, sq_J, weight, i);
			psi[i] += d1 * p[i];
		}
		k_times_psi(psi, phi, mass, nnt, ns, nn);
		apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
		for(i = 0; i < ns; i++){
			p[i] += c12 * p_dot(phi, tanh_Jphi, nnt, ns, nn, h, sq_J, weight, i);
		}
	}
	for(i = 0; i < ns; i++){
		p[i] -= c1 * p_dot(phi, tanh_Jphi, nnt, ns, nn, h, sq_J, weight, i);
	}
}

void leap_frog_parallel(double *psi, double *d_psi0, double *p, double weight, double h, double sq_J, double *tanh_averages, unsigned nmd, double dt, unsigned truncation_order, int up){
	const double psi0=*psi, dt_half=0.5*dt;
	double d_psi=*d_psi0, p0=*p;
	unsigned k;

	if(up == -1){
		p0 += p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order)*dt_half;
		printf("%g\t", p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order));
		d_psi += p0*dt_half;
	}
	for(k = 0; k < nmd; k++){
		d_psi += p0*dt_half;
		p0 += p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order)*dt;
		d_psi += p0*dt_half;
	}
	if(up == 1){
		d_psi += p0*dt_half;
		p0 += p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order)*dt_half;
		printf("%g\t", p_dot_parallel(psi0, d_psi, weight, h, sq_J, tanh_averages, truncation_order));
	}

	*d_psi0 = d_psi;
	*p = p0;
}

void leap_frog_split(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, unsigned sub_nmd, double sub_dt, unsigned truncation_order){
	const double dt_half=0.5*dt;
	const double weight = nn+mass, weight_parallel=2*nn+mass;
	double psi0, d_psi, d_psi_old, p0, shift;
	double force_sum;
	double *tanh_averages = malloc(truncation_order*sizeof(double));
	unsigned i, k;

	psi0 = sum_vector(psi, ns)/ns;
	p0 = sum_vector(p, ns)/ns;
	//expand_tanh(tanh_Jphi, tanh_averages, ns, truncation_order);
	for(i = 0; i < ns; i++) p[i] -= p0;

	for(k = 0; k < nmd; k++){
		//printf("%g\t%g\t%g\t%g\n", psi0, sum_vector(psi, ns)/ns, p0, sum_vector(p, ns)/ns);
		for(i = 0; i < ns; i++) psi[i] += p[i]*dt_half;
		update_fields_split(psi, phi, tanh_Jphi, tanh_averages, mass, sq_J, nnt, ns, nn, truncation_order);
		d_psi = 0;
		leap_frog_parallel(&psi0, &d_psi, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd/2, sub_dt, truncation_order, sub_nmd%2);
		for(i = 0; i < ns; i++) psi[i] += d_psi;
		d_psi_old = d_psi;
		leap_frog_parallel(&psi0, &d_psi, &p0, weight_parallel, h, sq_J, tanh_averages, sub_nmd/2, sub_dt, truncation_order, -1*(sub_nmd%2));
		psi0 += d_psi_old;
		d_psi -= d_psi_old;
		update_fields(psi, phi, tanh_Jphi, mass, sq_J, nnt, ns, nn);
		shift = weight_parallel*(-psi0 + sq_J*sum_vector(tanh_Jphi, ns)/ns);
		force_sum = 0;
		for(i = 0; i < ns; i++){
			p[i] += p_dot_orthogonal(phi, tanh_Jphi, nnt, ns, nn, sq_J, weight, i, shift)*dt;
			force_sum += pow(p_dot_orthogonal(phi, tanh_Jphi, nnt, ns, nn, sq_J, weight, i, shift),2);
			psi[i] += p[i]*dt_half + d_psi;
		}
		psi0 += d_psi;
		printf("%g\n", sqrt(force_sum/ns));
	}

	update_fields(psi, phi, tanh_Jphi, mass, sq_J, nnt, ns, nn);
	for(i = 0; i < ns; i++) p[i] += p0;
	free(tanh_averages);
}

void leap_frog(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt){
	const double weight = nn+mass;
	unsigned i, k;

	for(i = 0; i < ns; i++) psi[i] += dt/2*p[i];
	for(k = 0; k < nmd; k++){
		k_times_psi(psi, phi, mass, nnt, ns, nn);
		apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
		for(i = 0; i < ns; i++){
			p[i] += dt*p_dot(phi, tanh_Jphi, nnt, ns, nn, h, sq_J, weight, i);
			psi[i] += dt*p[i];
		}
	}
	for(i = 0; i < ns; i++){
		psi[i] -= dt/2*p[i];
	}
	k_times_psi(psi, phi, mass, nnt, ns, nn);
	apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
}

void integrate(double *psi, double *p, double *phi, double *tanh_Jphi, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, char *integrator){
	if(!strcmp(integrator, "leapfrog")){
		leap_frog(psi, p, phi, tanh_Jphi, nnt, ns, nn, h, sq_J, mass, nmd, dt);
	}else if(!strcmp(integrator, "4th")){
		evolve_Omelyan_4th(psi, p, phi, tanh_Jphi, nnt, ns, nn, h, sq_J, mass, nmd, dt);
	}else{
		char mode[100];
		unsigned sub_nmd, truncation_order;
		double sub_traj;
		int check = sscanf(integrator, "%s %u %lg %u ", mode, &sub_nmd, &sub_traj, &truncation_order);
	   	if(check == 4 && !strcmp(mode, "split-leapfrog")){
			leap_frog_split(psi, p, phi, tanh_Jphi, nnt, ns, nn, h, sq_J, mass, nmd, dt, sub_nmd, sub_traj*dt/sub_nmd, truncation_order);
		}else if(check == 4 && !strcmp(mode, "split-4th")){
			evolve_Omelyan_4th_split(psi, p, phi, tanh_Jphi, nnt, ns, nn, h, sq_J, mass, nmd, dt, sub_nmd, sub_traj*dt/sub_nmd, truncation_order);
		}else{
			printf("Integrator \"%s\" is not defined!\n", integrator);
			exit(0);
		}
	}
}

short trajectory(double *psi, double *p, unsigned *nnt, unsigned ns, unsigned nn, double h, double sq_J, double mass, unsigned nmd, double dt, double *psi_bar, char *integrator, gsl_rng *r){
	unsigned i;
	double psi_bar_old = *psi_bar;
	double energy_old, energy;
	double *phi = p+ns;
	double *tanh_Jphi = p+2*ns;
	double *psi_old = p+3*ns;
	double *phi_old = p+4*ns;
	double *tanh_Jphi_old = p+5*ns;

	k_times_psi(psi, phi, mass, nnt, ns, nn);
	apply_tanh_sq_J(phi, tanh_Jphi, sq_J, ns);
	for(i = 0; i < ns; i++){
		p[i] = gsl_ran_ugaussian(r);
		psi_old[i] = psi[i];
		phi_old[i] = phi[i];
		tanh_Jphi_old[i] = tanh_Jphi[i];
	}

	energy_old = hamilton(psi, phi, p, *psi_bar, sq_J, h, ns);
	integrate(psi, p, phi, tanh_Jphi, nnt, ns, nn, h, sq_J, mass, nmd, dt, integrator);
	*psi_bar = sum_vector(psi, ns)/ns;
	energy = hamilton(psi, phi, p, *psi_bar, sq_J, h, ns);
	//printf("%g\t%g\n", energy, energy_old);

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
	double deltaE = 2*ns*h/sq_J*(*psi_bar);
	if(exp(-deltaE) > gsl_rng_uniform(r)){
		for(unsigned i = 0; i < ns; i++) psi[i] *= -1;
		*psi_bar *= -1;
	}
}
