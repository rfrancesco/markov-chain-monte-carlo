#ifndef SCALAR1D_H
#define SCALAR1D_H
#include "../rng/ran2_class.hpp"

class Scalar1D {
  private:
    double ** lattice;
    unsigned int Nx;
    unsigned int Nt;
    double mass;
    double mass2;
    double m2p4;
    unsigned int * tnext;
    unsigned int * tprev;
    unsigned int * xnext;
    unsigned int * xprev;
    double Force(unsigned int, unsigned int) const;
    // Measurements buffer
    double m2phi2 = 0;
    double xdphi2 = 0;
    double tdphi2 = 0;
  public:
    explicit Scalar1D(unsigned int, unsigned int, double);
    ~Scalar1D();
    // Direct access to Lattice
    const double& operator() (unsigned int, unsigned int) const;
    double& operator() (unsigned int, unsigned int);
    // Measurements
    void CalculateObservables(); // Calculating observables separately is expensive
    double M2Phi2() const;
    double xdPhi2() const;
    double tdPhi2() const;
    double Energy();
    double TraceAnomaly();
    // Markov steps
    void HeatbathSweep();
    void OverrelaxationSweep();
    unsigned int HMC();
    // Debug
    void _PrintGeometry() const;
};


#endif
