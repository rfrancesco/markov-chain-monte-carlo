#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran2.c"


short random_spin() {
	double x = ran2();
	if (x > 0.5)
		return 1;
	else
		return -1;
}

short ** create_lattice(unsigned int l) {
	short ** matrix = malloc(l * sizeof(short*));
	unsigned int i;
	for (i = 0; i < l; i++)
		matrix[i] = malloc(l * sizeof(short));
	return matrix;
}

void initialize_lattice(short ** lattice, unsigned int l) {
	unsigned int i, j;
	for (i = 0; i < l; i++)
		for (j = 0; j < l; j++)
			lattice[i][j] = random_spin();
}

void printf_matrix(short ** lattice, unsigned int l) {
	unsigned int i, j;
	for (i = 0; i < l; i++) {
		for (j = 0; j < l; j++) 
			printf("%d\t", lattice[i][j]);
		printf("\n");
	}
}

int mod(int a, int b) {		// % non implementa il modulo ma il resto
	int r = a % b;
	return r < 0 ? r + b : r;
}

double energy_star(short ** lattice, unsigned int l, unsigned int i, unsigned int j) {
	double e = 0;
	e += lattice[mod(i+1, l)][j];
	e += lattice[mod(i-1, l)][j];
	e += lattice[i][mod(j+1, l)];
	e += lattice[i][mod(j-1, l)];
	e *= - lattice[i][j];
	return e;
}

void metropolis(short ** lattice, unsigned int l, double beta, unsigned int * accepted) {
	unsigned int i, j;
	double x;
	i = floor(ran2() * l);
	j = floor(ran2() * l);

	double p;
	p = exp(2 * beta * energy_star(lattice, l, i, j));
	x = ran2();
	if (x < p) {
		lattice[i][j] *= -1;
		(*accepted)++;
	}
}

double energy(short ** lattice, unsigned int l) {
	double e = 0;
	unsigned int i,j;
	for (i = 0; i < l; i++)
		for (j = 0; j < l; j++)
			e += energy_star(lattice, l, i, j);
	e /= 2;
	e /= (l*l);
	return e;
}

double magnetization(short ** lattice, unsigned int l) {
	double m = 0;
	unsigned int i,j;
	for (i = 0; i < l; i++)
		for(j = 0; j < l; j++)
			m += lattice[i][j];
	m /= (l*l);
	return m;
}

int main() {

	ran2_init();
	double beta = 0.3;
	unsigned int l = 10;
	unsigned int n_measures = 100000;
	unsigned int n_skip = 10;
	n_skip *= l * l;
	short ** lattice = create_lattice(l);
	initialize_lattice(lattice, l);

	//printf_matrix(lattice, l);	
	
	unsigned int i, j;
	unsigned int accept = 0; 	// accettanza (da normalizzare)


	printf("#n\tEnergy\tMagnetization\n");
	for (i = 0; i < n_measures; i++) {
		for (j = 0; j < n_skip; j++)
			metropolis(lattice, l, beta, &accept);
		printf("%u\t%f\t%f\n", i, energy(lattice, l), magnetization(lattice, l));
	}
	printf("#Acceptance: %u out of %u (%f)\n", accept, n_measures*n_skip, ((float) accept) / (n_measures*n_skip));
	ran2_save();
	return 0;
}
