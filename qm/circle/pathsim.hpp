#ifndef PATHSIM_H
#define PATHSIM_H
#include <cmath>
#include "s1.hpp"
#include "../../rng/ran.hpp"

class PathSim {
    private:
        S1 * p;
        unsigned int size;
        unsigned int* prev;
        unsigned int* next;
    public:
        // Constructor and destructor
        explicit PathSim(unsigned int);
        ~PathSim();
        // Direct access to path
        const S1& operator[] (unsigned int) const;
        S1& operator[] (unsigned int);
        // Index manipulation
        unsigned int p1(unsigned int) const;
        unsigned int m1(unsigned int) const;
        unsigned int add(unsigned int, unsigned int) const;
        unsigned int sub(unsigned int, unsigned int) const;
        //void print(double) const;
        // Observables
        double winding() const;
        double slab_winding(unsigned int) const;
        double cos_avg() const;
        double cos_correlator(unsigned int) const;
        double cos2_avg() const;
        double cos2_correlator(unsigned int) const;
        // Markov steps - Choose your weapon!
        int MetropolisSweep(double, double);
        int TailorStep(double);
};

#endif
