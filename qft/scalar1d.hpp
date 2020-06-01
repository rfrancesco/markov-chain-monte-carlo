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
    // Periodic boundary conditions
    unsigned int * tnext;
    unsigned int * tprev;
    unsigned int * xnext;
    unsigned int * xprev;
    // Force() needed by Markov steps functions
    double Force(unsigned int, unsigned int) const;
    // Measurements buffer
    double m2phi2 = 0;
    double xdphi2 = 0;
    double tdphi2 = 0;
  public:
    // Constructor/Destructor
    explicit Scalar1D(unsigned int, unsigned int, double);
    ~Scalar1D();
    // Direct access to Lattice
    const double& operator() (unsigned int, unsigned int) const;
    double& operator() (unsigned int, unsigned int);
    // Measurements
    void CalculateObservables(); // Calculating these three observables separately is expensive
    double M2Phi2() const;
    double xdPhi2() const;
    double tdPhi2() const;
    // Correlators
    double FTCorrelator(unsigned int, double, double) const;
    // Markov steps
    void HeatbathSweep();
    void OverrelaxationSweep();
    // Debug
    void _PrintGeometry() const;
};


#endif
