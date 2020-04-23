#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../rng/ran2.c"

typedef struct path_s {
	double * p;
	unsigned int len;
} path_t;

// Funzioni di inizializzazione del cammino
path_t * create_path(unsigned int n) {			
	path_t * path = malloc(sizeof(path_t));
	if (!path)
		return NULL;
	path->p = malloc(n * sizeof(double));	
	if (!path->p) {
		free(path);
		return NULL;
	}
	path->len = n;
	return path;
}

void free_path(path_t * path) {
	free(path->p);
	free(path);
}

void initialize_path(path_t * path) {
	// per ora mi limito a inizializzare a zero
	for (unsigned int i = 0; i < path->len; i++)
		path->p[i] = 0;
}

// Funzione di debug sul cammino
void print_path(path_t * path) {
	for (unsigned int i = 0; i < path->len; i++)
		printf("%f\t", path->p[i]);
	printf("\n");
}

// In C, % non implementa il modulo => workaround
// Da sostituire con lookup
unsigned int mod(int a, int b) {
	int r = a % b;
	return r < 0 ? r + b : r;
}

// Correlatore <x(0)x(τ)>
double x_corr(path_t * path, unsigned int k) {
	double c = 0;
	for (unsigned int i = 0; i < path->len; i++)
		c += path->p[i] * path->p[mod(i+k, path->len)];
	c /= path->len;
	return c;
}


void metropolis(path_t * path, double eta, double delta, unsigned int *accepted) {
	for (unsigned int i = 0; i < path->len; i++) { 			//spazza il cammino, rispetta il bilancio dettagliato
		double yp = path->p[i] + delta*(2*ran2() - 1);		// y[i] di prova

		double r = exp( - (pow(yp,2) - pow(path->p[i],2))*((eta/2) + (1 / eta)) + (1/eta) * (yp - path->p[i])*(path->p[mod(i+1, path->len)] + path->p[mod(i-1, path->len)]));

		if (r > 1) {						// Metropolis test
			(*accepted)++;
			path->p[i] = yp;
		}
		else {
			double x = ran2();
			if (x < r) {
				(*accepted)++;
				path->p[i] = yp;
			}
		}	
	}
}

double x2(path_t * path) {
	double m = 0;
	for (unsigned int i = 0; i < path->len; i++)
		m += pow(path->p[i],2);
	m /= path->len;
	return m;
}

double dx2(path_t * path) {
	double m = 0;
	for (unsigned int i = 0; i < path->len; i++)
		m += pow((path->p[i] - path->p[mod(i-1, path->len)]), 2);
	m /= path->len;
	return m;
}


int main(int argc, char *argv[]) {
	if (argc < 4) {
		printf("Syntax is %s [N*eta] [N] [n_measures] [n_skip] \n", argv[0]);
	    	return EXIT_FAILURE;
	}	
	ran2_init();
	//double eta = 0.03;
	double Neta = atof(argv[1]);
	unsigned int N = atoi(argv[2]);
	double eta = Neta / N;
	printf("#Beginning simulation with Neta = %f, N = %u, eta = %f\n", Neta, N, eta);
	unsigned int n_measures = atoi(argv[3]);
	unsigned int n_skip = atoi(argv[4]);
	unsigned int accept = 0;
	double delta = sqrt(eta);
	//double delta = atof(argv[3]);
	path_t * path = create_path(N);
	initialize_path(path);

	printf("#y^2\tdy^2\t<x(τ)x(0)>\n");
	for (unsigned int i = 0; i < n_measures; i++) {
		for (unsigned int j = 0; j < n_skip; j++)
			metropolis(path, eta, delta, &accept);
		printf("%f\t%f\t", x2(path), dx2(path));
		for (unsigned int t = 0; t < path->len; t++)
			printf("%f\t", x_corr(path, t));
		printf("\n");
	}

	printf("#Acceptance: %u out of %u (%f)\n", accept, n_measures*n_skip*N, ((float) accept) / (n_measures*n_skip*N));
	free_path(path);
	ran2_save();
	return EXIT_SUCCESS;
}
