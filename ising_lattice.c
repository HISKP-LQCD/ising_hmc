#include "ising_lattice.h"

unsigned *construct_lattice(char *geometry, FILE *input, unsigned *ns, unsigned *nn, double **s){
	unsigned *nnt;
	if(!strcmp(geometry, "cubic")){
		nnt = construct_cubic(input, ns, nn, s);
	}else if(!strcmp(geometry, "rectangular")){
		nnt = construct_rectangular(input, ns, nn, s);
	}else if(!strcmp(geometry, "triangular")){
		nnt = construct_triangular(input, ns, nn, s);
	}else{
		printf("The geometry \"%s\" is not known!\n", geometry);
		exit(0);
	}
	fclose(input);
	return nnt;
}

unsigned *construct_cubic(FILE *input, unsigned *ns, unsigned *nn, double **s){
	unsigned dim, l, d, i;
	unsigned *N, *nnt;
	char boundaries[100];
	int id;

	fscanf(input, "dimension=%u\n", &dim);
	fscanf(input, "length=%u\n", &l);
	fscanf(input, "boundaries=%s\n", boundaries);

	*nn = 2*dim;
	N = malloc(dim * sizeof(unsigned));
	N[0] = 1;
	for(d = 1; d < dim; d++){
		N[d] = N[d-1]*l;
	}
	*ns = N[d-1]*l;

	const unsigned size = *ns;
	const unsigned neighbours = *nn;
	*s = malloc((size+1) * sizeof(double));
	nnt = malloc(neighbours * size * sizeof(unsigned));

	if(!strcmp(boundaries, "periodic")){
		for(i = 0; i < size; i++){
			for(d = 0; d < dim; d++){
				id = (i/N[d]) % l;
				nnt[i*neighbours + 2*d] = i + ((id+1)%l - id)*N[d];
				nnt[i*neighbours + 2*d + 1] = i + ((id+l-1)%l - id)*N[d];
			}
		}
	}else if(!strcmp(boundaries, "open") || !strcmp(boundaries, "closed")){
		if(!strcmp(boundaries, "open")) (*s)[size] = 0;
		else (*s)[size] = 1;
		for(i = 0; i < size; i++){
			for(d = 0; d < dim; d++){
				id = (i/N[d]) % l;
				if(id < l-1) nnt[i*neighbours + 2*d] = i + N[d];
				else nnt[i*neighbours + 2*d] = size;
				if(id > 0) nnt[i*neighbours + 2*d + 1] = i - N[d];
				else nnt[i*neighbours + 2*d + 1] = size;
			}
		}
	}else{
		printf("Boundary conditions \"%s\" not known!\n", boundaries);
		exit(0);
	}

	free(N);
	return nnt;
}

unsigned *construct_rectangular(FILE *input, unsigned *ns, unsigned *nn, double **s){
	unsigned dim, d, i;
	unsigned *l, *N, *nnt;
	char boundaries[100];
	int id;

	fscanf(input, "dimension=%u\n", &dim);
	l = malloc(dim * sizeof(unsigned));
	N = malloc(dim * sizeof(unsigned));
	fscanf(input, "lengths=%u ", l);
	for(d = 1; d < dim; d++){
		fscanf(input, "%u ", l+d);
	}
	fscanf(input, "boundaries=%s\n", boundaries);

	*nn = 2*dim;
	N[0] = 1;
	for(d = 1; d < dim; d++){
		N[d] = N[d-1]*l[d-1];
	}
	*ns = N[d-1]*l[d-1];

	const unsigned size = *ns;
	const unsigned neighbours = *nn;
	*s = malloc((size+1) * sizeof(double));
	nnt = malloc(neighbours * size * sizeof(unsigned));

	if(!strcmp(boundaries, "periodic")){
		for(i = 0; i < size; i++){
			for(d = 0; d < dim; d++){
				id = (i/N[d]) % l[d];
				nnt[i*neighbours + 2*d] = i + ((id+1)%l[d] - id)*N[d];
				nnt[i*neighbours + 2*d + 1] = i + ((id+l[d]-1)%l[d] - id)*N[d];
			}
		}
	}else if(!strcmp(boundaries, "open") || !strcmp(boundaries, "closed")){
		if(!strcmp(boundaries, "open")) (*s)[size] = 0;
		else (*s)[size] = 1;
		for(i = 0; i < size; i++){
			for(d = 0; d < dim; d++){
				id = (i/N[d]) % l[d];
				if(id < l[d]-1) nnt[i*neighbours + 2*d] = i + N[d];
				else nnt[i*neighbours + 2*d] = size;
				if(id > 0) nnt[i*neighbours + 2*d + 1] = i - N[d];
				else nnt[i*neighbours + 2*d + 1] = size;
			}
		}
	}else{
		printf("Boundary conditions \"%s\" not known!\n", boundaries);
		exit(0);
	}

	free(l);
	free(N);
	return nnt;
}

unsigned *construct_triangular(FILE *input, unsigned *ns, unsigned *nn, double **s){
	// constructs a parallelogram with periodic boundaries
	unsigned l1, l2, i;
	unsigned *nnt;
	int i1, i2;

	fscanf(input, "lengths=%u %u\n", &l1, &l2);

	*nn = 6;
	*ns = l1*l2;

	const unsigned size = *ns;
	const unsigned neighbours = *nn;
	*s = malloc(size * sizeof(double));
	nnt = malloc(neighbours * size * sizeof(unsigned));

	for(i = 0; i < size; i++){
		i1 = i%l1;
		i2 = i/l1;
		nnt[i*neighbours] = i + (i1+1)%l1 - i1; // right
		nnt[i*neighbours+1] = i + (i1+l1-1)%l1 - i1; // left
		nnt[i*neighbours+2] = i + ((i2+1)%l2 - i2)*l1; // top right
		nnt[i*neighbours+3] = i + ((i2+l2-1)%l2 - i2)*l1; // bottom left
		nnt[i*neighbours+4] = nnt[i*neighbours+2] + (i1+l1-1)%l1 - i1; // top left
		nnt[i*neighbours+5] = nnt[i*neighbours+3] + (i1+1)%l1 - i1; // bottom right
	}

	return nnt;
}