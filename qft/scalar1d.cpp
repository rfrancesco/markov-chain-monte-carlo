#include "scalar1d.hpp"
#include <cmath>

extern Ran2 ran;

Scalar1D::Scalar1D(unsigned int Nx, unsigned int Nt, double mass) {
    this->Nx = Nx;
    this->Nt = Nt;
    this->mass = mass;
    mass2 = pow(mass,2);
    m2p4 = mass2 + 4;
    // Ignoring memory issues for now
    lattice = new double*[Nx];
    for (unsigned int x = 0; x < Nx; ++x)
        lattice[x] = new double[Nt];
    // Setting up periodic conditions...
    xnext = new unsigned int[Nx];
    xprev = new unsigned int[Nx];
    tnext = new unsigned int[Nt];
    tprev = new unsigned int[Nt];
    for (unsigned int x = 0; x < Nx; ++x) {
        xnext[x] = x + 1;
        xprev[x] = x - 1;
    }
    xnext[Nx - 1] = 0;
    xprev[0] = Nx - 1;
    for (unsigned int t = 0; t < Nt; ++t) {
        tnext[t] = t + 1;
        tprev[t] = t - 1;
    }
    tnext[Nt - 1] = 0;
    tprev[0] = Nt - 1;
    // Initialization to zero
    for (unsigned int x = 0; x < Nx; ++x)
        for (unsigned int t = 0; t < Nt; ++t)
            lattice[x][t] = 0;
}

Scalar1D::~Scalar1D() {
    // Cleanup
    delete [] xnext;
    delete [] xprev;
    delete [] tnext;
    delete [] tprev;
    for (unsigned int x = 0; x < Nx; ++x)
        delete [] lattice[x];
    delete [] lattice;
}


// Direct access to Lattice
const double& Scalar1D::operator()(unsigned int x, unsigned int t) const {
    return lattice[x][t];
}

double& Scalar1D::operator()(unsigned int x, unsigned int t) {
    return lattice[x][t];
}

// Measurements

double Scalar1D::TraceAnomaly() {
    // returns (ε - p)/T²
    double phi2_avg = 0;
    for (unsigned int x = 0; x < Nx; ++x)
        for (unsigned int t = 0; t < Nt; ++t)
            phi2_avg += pow(lattice[x][t], 2);
    return phi2_avg * Nx * Nt / mass2;
}

// Markov steps

double Scalar1D::Force(unsigned int x, unsigned int t) const {
    double f = lattice[x][tnext[t]] + lattice[x][tprev[t]];
    f += lattice[xnext[x]][t] + lattice[xprev[x]][t];

    return f;
}

void Scalar1D::HeatbathSweep() {
    for (unsigned int x = 0; x < Nx; ++x)
        for (unsigned int t = 0; t < Nt; ++t) {
            double mean = Force(x,t)/m2p4;
            double variance = 1/m2p4;

            lattice[x][t] = ran.gaussian(mean, variance);
        }
}

void Scalar1D::OverrelaxationSweep() {
    for (unsigned int x = 0; x < Nx; ++x)
        for (unsigned int t = 0; t < Nt; ++t) {
            double mean = Force(x,t)/m2p4;
            lattice[x][t] = (2 * mean) - lattice[x][t];
        }
}
