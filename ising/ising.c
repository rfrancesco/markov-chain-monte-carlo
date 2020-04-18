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
    unsigned int * next;
    unsigned int * prev;
    unsigned int size;
} lattice_t;

void ungracefully_handle_memory_issues() {
    printf("Fatal error: Memory allocation failure.\n");
    fprintf(stderr, "Fatal error: Memory allocation failure.\n");
    exit(EXIT_FAILURE);
}

short ** create_matrix(unsigned int l) {
    short ** matrix = malloc(l * sizeof(short*));
    if (!matrix)
        ungracefully_handle_memory_issues();
    unsigned int i;
    for (i = 0; i < l; i++) {
        matrix[i] = malloc(l * sizeof(short));
        if (!matrix[i])
            ungracefully_handle_memory_issues();
    }
    return matrix;
}

void free_matrix(short ** matrix, unsigned int l) {
    for (unsigned int i = 0; i < l; i++)
        free(matrix[i]);
    free(matrix);
}

lattice_t * create_lattice(unsigned int l) {
    lattice_t * lat = malloc(sizeof(lattice_t));
    lat->matrix = create_matrix(l);
    lat->size = l;
    lat->next = malloc(l*sizeof(unsigned int));
    lat->prev = malloc(l*sizeof(unsigned int));
    if (!lat->prev || !lat->next)
        ungracefully_handle_memory_issues();
    for (unsigned int i = 0; i < lat->size; i++) {
        lat->next[i] = i+1;
        lat->prev[i] = i-1;
    }
    lat->next[l-1] = 0;
    lat->prev[0] = l-1;
    return lat;
}

void free_lattice(lattice_t * lat) {
    free(lat->next);
    free(lat->prev);
    free_matrix(lat->matrix, lat->size);
    free(lat);
}

void initialize_lattice(lattice_t * lat, unsigned int init) {
    unsigned int i, j;
    switch (init) {
        case 0:	for (i = 0; i < lat->size; i++)
                    for (j = 0; j < lat->size; j++)
                        lat->matrix[i][j] = 1;
                printf("#Lattice initialized with 1\n");
                break;
        case 1:	for (i = 0; i < lat->size; i++)
                    for (j = 0; j < lat->size; j++)
                        lat->matrix[i][j] = random_spin();
                printf("#Lattice initialized with random spins\n");
                break;
        default:    printf("Fatal error: Illegal value for init\n");
                    fprintf(stderr, "Fatal error: Illegal value for init\n");
                    exit(EXIT_FAILURE);
    }
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

double force(lattice_t * lat, unsigned int i, unsigned int j) {
    // Using mod() is convenient but wastes time (~15%)
    // A lookup table would be a more efficient approach (see ising.cpp)
    double e = 0;
    e += lat->matrix[lat->next[i]][j];
    e += lat->matrix[i][lat->next[j]];
    e += lat->matrix[lat->prev[i]][j];
    e += lat->matrix[i][lat->prev[j]];
    return e;
}

void metropolis(lattice_t * lat, double beta, double extfield, unsigned long int * accepted) {
    unsigned int i, j;
    i = floor(ran2() * lat->size);
    j = floor(ran2() * lat->size);

    double p = exp(- 2 * beta * lat->matrix[i][j] * (force(lat, i, j) + extfield));
    if (p > 1) {
        lat->matrix[i][j] *= -1;
        (*accepted)++;
    } 
    else {
        double x = ran2();
        if (x < p) {
            lat->matrix[i][j] *= -1;
            (*accepted)++;
        }
    }
}

double energy(lattice_t * lat) {
    double e = 0;
    unsigned int i,j;
    for (i = 0; i < lat->size; i++)
        for (j = 0; j < lat->size; j++)
            e -= lat->matrix[i][j]*force(lat, i, j);
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

int main(int argc, char** argv) {

    if (argc != 7) {
        printf("Usage: n_measures, n_skip, L, beta, extfield, init (0: cold, 1: hot)\n");
        return EXIT_FAILURE;
    }
        
    ran2_init();
    unsigned int n_measures = atoi(argv[1]);
    unsigned int n_skip = atoi(argv[2]);
    unsigned int l = atoi(argv[3]);
    double beta = atof(argv[4]);
    double extfield = atof(argv[5]);
    unsigned int init = atoi(argv[6]);
    printf("#Beginning simulation: L = %u, beta = %f, extfield = %f\n", l, beta, extfield);
    printf("#%u measures every %u sweeps\n", n_measures, n_skip);
    
    // n_skip is measured in sweeps
    n_skip *= l * l;
    lattice_t * lattice = create_lattice(l);
    initialize_lattice(lattice, init);

    //printf_matrix(lattice, l);	
    
    unsigned int i, j;
    unsigned long int accepted = 0; 	// counts accepted steps 


    printf("#Energy\tMagnetization\n");
    for (i = 0; i < n_measures; i++) {
        for (j = 0; j < n_skip; j++)
            metropolis(lattice, beta, extfield, &accepted);
        printf("%f\t%f\n", energy(lattice), magnetization(lattice));
    }
    printf("#Acceptance: %u out of %u (%f)\n", accepted, n_measures*n_skip, ((float) accepted) / (n_measures*n_skip));
    ran2_save();
    free_lattice(lattice);
    return 0;
}
