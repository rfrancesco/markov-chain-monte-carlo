#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../rng/ran2.c"

double gauss(double x, double m, double s2) {
  // this is not normalized, but it doesn't matter
  return exp(-(x-m)*(x-m)/(2*s2));
}

double metropolis(double x, double delta, double m, double s2, unsigned int * accepted) {
  double xp = x + delta * (2*ran2() - 1);
  double r = gauss(xp,m,s2) / gauss(x,m,s2);
  if (r > 1) {
   (*accepted)++;
   return xp;
  }
  else {
    double p = ran2();
    if (p < r) {
      (*accepted)++;
      return xp;
    }
    else
      return x;
  }
}


int main(int argc, char** argv ) {

  if (argc < 3) {
    printf("Usage: %s [average] [variance]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  ran2_init();
  // Read mean and average from argv
  double m = atof(argv[1]);
  double s2 = atof(argv[2]);
  double delta = s2/10;
  unsigned long int Nmax = 110000;
  unsigned int accepted = 0;
  printf("#Simulating m = %f, s2 = %f\n", m, s2);
  double x = m - 2*s2;

  for (unsigned long int i = 0; i < Nmax; i++) {
    x = metropolis(x, delta, m, s2, &accepted);
    printf("%f\n", x);
  }

  printf("#Acceptance: %f\n", ((float) accepted) / Nmax); 

  ran2_save();
  return EXIT_SUCCESS;
}
