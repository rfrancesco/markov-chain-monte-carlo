#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <cmath>
#include <random>

#include "s1.hpp"
#include "../../rng/ran2_class.hpp"
#include "pathsim.hpp"

const double PI2 = pow(M_PI, 2);
std::random_device urandom("/dev/urandom");
Ran2 ran(urandom());

int main(int argc, char** argv) {

    if (argc != 5) {
        cout << "Usage: 4 arguments: Neta, N, n_measures, n_skip" << endl;
        return EXIT_FAILURE;
    }

    const double Neta = atof(argv[1]); // è essenzialmente la temperatura: Nη = β/mL^2
    const unsigned int N = atoi(argv[2]);
    const double eta = Neta/N;
    const unsigned long int n_measures = atoi(argv[3]);
    const unsigned long int n_skip = atoi(argv[4]);

    const double delta = 0.5;

    unsigned long int accepted = 0;

    cout << "#Simulation with Νη = β/mL^2 = " << Neta << ", η = " << eta << ", N = " << N << endl;
    cout << "#Taking " << n_measures << " measures every " << n_skip << "sweeps" << endl;
    cout << "#Local Metropolis algorithm with delta = " << delta << endl;

    PathSim* p = new PathSim(N);

    for (unsigned int i = 0; i < n_measures; ++i) {
        // Metropolis steps
        for (unsigned int j = 0; j < n_skip; ++j)
            accepted += p->MetropolisSweep(eta, delta);
        // Print measures
        printf("%.3f\n", p->winding());
    }

    cout << "#Acceptance: " << accepted << " out of " << N*n_measures*n_skip << " (" << static_cast<double>(accepted)/(N*n_measures*n_skip) << ")" << endl;

    return 0;
}
