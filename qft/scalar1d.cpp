#include "scalar1d.hpp"
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;
#include <complex>
using std::complex;

extern Ran2 ran;

Scalar1D::Scalar1D(unsigned int Nx, unsigned int Nt, double mass) {
    this->Nx = Nx;
    this->Nt = Nt;
    this->mass = mass;
    mass2 = mass * mass;
    m2p4 = mass2 + 4;
    // Ignoring memory issues for now
    lattice = new double*[Nx];
    for (unsigned int x = 0; x < Nx; ++x)
        lattice[x] = new double[Nt];
    // Setting up periodic conditions...
    // Space
    xnext = new unsigned int[Nx];
    xprev = new unsigned int[Nx];
    for (unsigned int x = 0; x < Nx; ++x) {
        xnext[x] = x + 1;
        xprev[x] = x - 1;
    }
    xnext[Nx - 1] = 0;
    xprev[0] = Nx - 1;
    // Time
    tnext = new unsigned int[Nt];
    tprev = new unsigned int[Nt];
    for (unsigned int t = 0; t < Nt; ++t) {
        tnext[t] = t + 1;
        tprev[t] = t - 1;
    }
    tnext[Nt - 1] = 0;
    tprev[0] = Nt - 1;
    // Initialization to zero
    for (unsigned int x = 0; x < Nx; ++x)
        for (unsigned int t = 0; t < Nt; ++t)
            lattice[x][t] = ran.doub();
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
double Scalar1D::M2Phi2() const {
    return m2phi2;
}

double Scalar1D::xdPhi2() const {
    return xdphi2;
}

double Scalar1D::tdPhi2() const {
    return tdphi2;
}

void Scalar1D::CalculateObservables() {
    m2phi2 = 0;
    xdphi2 = 0;
    tdphi2 = 0;
    for (unsigned int x = 0; x < Nx; ++x)
        for (unsigned int t = 0; t < Nt; ++t) {
            m2phi2 += pow(lattice[x][t],2);
            tdphi2 += pow(lattice[x][tnext[t]] - lattice[x][t], 2);
            xdphi2 += pow(lattice[xnext[x]][t] - lattice[x][t], 2);
        }
    m2phi2 *= mass2 / (Nx * Nt);
    tdphi2 /= Nx * Nt;
    xdphi2 /= Nx * Nt;
}

// Correlators
double Scalar1D::FTCorrelator(unsigned int k, double p, double q) const {
    complex<double> result = 0;
    unsigned int t0 = 0;
    unsigned int tk = k;
    for (t0 = 0; t0 < Nt; ++t0) {
        complex<double> ftx = 0;
        complex<double> fty = 0;
        for (unsigned int x = 0; x < Nx; ++x) {
            ftx += std::polar(1., -q*x) * lattice[x][tk];
            fty += std::polar(1., p*x) * lattice[x][t0];
        }
        result += ftx * fty / static_cast<double>(Nx);
        tk = tnext[tk];
    }
    result /= Nt;
    return abs(result);
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

// Debug
void Scalar1D::_PrintGeometry() const {
    cout << "X" << endl;
    for (unsigned int x = 0; x < Nx; ++x)
        cout << xprev[x] << x << xnext[x] << endl;

    cout << "T" << endl;
    for (unsigned int t = 0; t < Nt; ++t)
        cout << tprev[t] << t << tnext[t] << endl;
}
