#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../rng/ran2.c"


short random_spin() {
	double x = ran2();
	if (x > 0.5)
		return 1;
	else
		return -1;
}

typedef struct lattice_s {
	short ** matrix;
	unsigned int size;
} lattice_t;


short ** create_matrix(unsigned int l) {
	short ** matrix = malloc(l * sizeof(short*));
	unsigned int i;
	for (i = 0; i < l; i++)
		matrix[i] = malloc(l * sizeof(short));
	return matrix;
}

lattice_t * create_lattice(unsigned int l) {
	lattice_t * lat = malloc(sizeof(lattice_t));
	lat->matrix = create_matrix(l);
	if (lat->matrix == NULL) {
		free(lat);
		return NULL;
	}
	else {
		lat->size = l;
		return lat;
	}
}

void initialize_lattice(lattice_t * lat) {
	unsigned int i, j;
	for (i = 0; i < lat->size; i++)
		for (j = 0; j < lat->size; j++)
			lat->matrix[i][j] = random_spin();
}

/*void printf_matrix(short ** lattice, unsigned int l) {
	unsigned int i, j;
	for (i = 0; i < l; i++) {
		for (j = 0; j < l; j++) 
			printf("%d\t", lattice[i][j]);
		printf("\n");
	}
}*/

int mod(int a, int b) {		// C % is remainder, not modulus
	int r = a % b;
	return r < 0 ? r + b : r;
}

double energy_neighbors(lattice_t * lat, unsigned int i, unsigned int j) {
	// Using mod() is convenient but wastes time (~15%)
	// A lookup table would be a more efficient approach (see ising.cpp)
	double e = 0;
	e += lat->matrix[mod(i+1, lat->size)][j];
	e += lat->matrix[mod(i-1, lat->size)][j];
	e += lat->matrix[i][mod(j+1, lat->size)];
	e += lat->matrix[i][mod(j-1, lat->size)];
	e *= - lat->matrix[i][j];
	return e;
}

void metropolis(lattice_t * lat, double beta, unsigned int * accepted) {
	unsigned int i, j;
	double x;
	i = floor(ran2() * lat->size);
	j = floor(ran2() * lat->size);

	double p;
	p = exp(2 * beta * energy_neighbors(lat, i, j));
	x = ran2();
	if (x < p) {
		lat->matrix[i][j] *= -1;
		(*accepted)++;
	}
}

double energy(lattice_t * lat) {
	double e = 0;
	unsigned int i,j;
	for (i = 0; i < lat->size; i++)
		for (j = 0; j < lat->size; j++)
			e += energy_neighbors(lat, i, j);
	e /= 2;
	e /= pow(lat->size,2);
	return e;
}

double magnetization(lattice_t * lat) {
	double m = 0;
	unsigned int i,j;
	for (i = 0; i < lat->size; i++)
		for(j = 0; j < lat->size; j++)
			m += lat->matrix[i][j];
	m /= pow(lat->size,2);
	return m;
}

int main() {

	ran2_init();
	double beta = 0.3;
	unsigned int l = 10;
	unsigned int n_measures = 100000;
	unsigned int n_skip = 100;

	n_skip *= l * l;
	lattice_t * lattice = create_lattice(l);
	initialize_lattice(lattice);

	//printf_matrix(lattice, l);	
	
	unsigned int i, j;
	unsigned int accept = 0; 	// accettanza (da normalizzare)


	printf("#Energy\tMagnetization\n");
	for (i = 0; i < n_measures; i++) {
		for (j = 0; j < n_skip; j++)
			metropolis(lattice, beta, &accept);
		printf("%f\t%f\n", energy(lattice), magnetization(lattice));
	}
	printf("#Acceptance: %u out of %u (%f)\n", accept, n_measures*n_skip, ((float) accept) / (n_measures*n_skip));
	ran2_save();
	return 0;
}
