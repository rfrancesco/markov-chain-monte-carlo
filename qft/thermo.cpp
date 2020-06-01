#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <cmath>
#include <random>

#include "../rng/ran2_class.hpp"
#include "scalar1d.hpp"

std::random_device urandom("/dev/urandom");
Ran2 ran(urandom());

int main(int argc, char** argv) {

    if (argc != 5) {
        cout << "Usage: 4 arguments: Nx, Nt, m, n_measures" << endl;
        return EXIT_FAILURE;
    }
    const unsigned int Nx = atoi(argv[1]);
    const unsigned int Nt = atoi(argv[2]);
    const double m = atof(argv[3]);
    const unsigned long int n_measures = atoi(argv[4]);

    cout << "#Simulation of (Nx, Nt) = " << Nx << ", " << Nt << ") with mass" << m << endl;
    cout << "#Taking " << n_measures << "  measures every HB+4OR" << endl;

    Scalar1D* field = new Scalar1D(Nx, Nt, m);

    for (unsigned int i = 0; i < n_measures; ++i) {
        field->HeatbathSweep();
        field->OverrelaxationSweep();
        field->OverrelaxationSweep();
        field->OverrelaxationSweep();
        field->OverrelaxationSweep();
        field->CalculateObservables();
        cout << field->M2Phi2() << "\t" << field->xdPhi2() << "\t" << field->tdPhi2() << endl;
    }

    cout << "#Simulation finished" << endl;
    return 0;

}
