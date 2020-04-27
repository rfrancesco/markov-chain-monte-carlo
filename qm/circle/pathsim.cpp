#include "pathsim.hpp"

// RNG defined as a global variable
// This is not ideal but ok for the purposes of my programs
// To be more general I'd probably have to define a Ran*
// and pass it to the constructor

extern Ran2 ran;

// Constructor and destructor
PathSim::PathSim(unsigned int s) {
    p = new S1[s];
    // Allocating and setting prev, next implementing boundary conditions 
    prev = new unsigned int[s];
    next = new unsigned int[s];
    for (unsigned int i = 0; i < s; i++) {
        next[i] = i + 1;
        prev[i] = i - 1;
    }
    next[s-1] = 0;
    prev[0] = s - 1;
    size = s;
    //for (unsigned int i = 0; i < size; i++)
    //    p[i] = ran.doub();
}

PathSim::~PathSim() {
    delete [] prev;
    delete [] next;
    delete [] p;
}

// Direct access to path
const S1& PathSim::operator[](unsigned int i) const {
    return p[i];
}

S1& PathSim::operator[](unsigned int i) {
    return p[i];
}

// Index manipulation
unsigned int PathSim::p1(unsigned int i) const {
    return next[i];
}

unsigned int PathSim::m1(unsigned int i) const {
    return prev[i];
}

unsigned int PathSim::add(unsigned int i, unsigned int k) const {
    while (k > 0) {
        i = p1(i);
        --k;
    }
    return i;
}

unsigned int PathSim::sub(unsigned int i, unsigned int k) const {
    while (k > 0) {
        i = m1(i);
        --k;
    }
    return i;
}

// Observables
double PathSim::winding() const {
    double q = 0;
    for (unsigned int i = 0; i < size; ++i)
        q += p[p1(i)] - p[i];
    return q;
}

double PathSim::slab_winding(unsigned int stop) const {
    double q = 0;
    for (unsigned int i = 0; i < stop; ++i)
        q += p[p1(i)] - p[i];
    return q;
}

double PathSim::cos_avg() const {
    double result = 0;
    // Calculate <x^2> as average over entire path
    for (unsigned int i = 0; i < size; ++i)
        result += cos(2*M_PI*p[i]);
    result /= size;
    return result;
}

double PathSim::cos_correlator(unsigned int k) const {
    double result = 0;
    unsigned int i = 0;
    unsigned int ik = add(i,k);
    // Calculate C(k) as average of C(0, k), C(1, 1+k) ... C(s-1, s-1+k)
    for (i = 0; i < size; ++i) {
        result += cos(2*M_PI*p[ik])*cos(2*M_PI*p[i]);
        ik = p1(ik);
    }
    result /= size;
    return result;
}

double PathSim::cos2_avg() const {
    double result = 0;
    // Calculate <x^2> as average over entire path
    for (unsigned int i = 0; i < size; ++i)
        result += cos(4*M_PI*p[i]);
    result /= size;
    return result;
}

double PathSim::cos2_correlator(unsigned int k) const {
    double result = 0;
    unsigned int i = 0;
    unsigned int ik = add(i,k);
    // Calculate C(k) as average of C(0, k), C(1, 1+k) ... C(s-1, s-1+k)
    for (i = 0; i < size; ++i) {
        result += cos(4*M_PI*p[ik])*cos(4*M_PI*p[i]);
        ik = p1(ik);
    }
    result /= size;
    return result;
}




/*  void PathSim::print(double Neta) const {
    printf("%.3f\t", winding());
    for (unsigned int i = 0; i < size/(7*Neta); i++)
        cout << cos_correlator(i) << "\t"; // it is already a connected correlator
    double cos2_avg2 = pow(cos2_avg(),2);
    for (unsigned int i = 0; i < size/(7*Neta); i++)
        cout << cos2_correlator(i) - cos2_avg2 << "\t"; // connected correlator
    cout << '\n';
} */



// Markov steps
int PathSim::MetropolisSweep(double eta, double delta) {
    unsigned int accepted = 0;
    for (unsigned int i = 0; i < size; i++) {
        double x = delta*(2*ran.doub() - 1);
        S1 yp = p[i] + x;

        double dS = pow( p[p1(i)] - yp , 2);
        dS += pow(yp - p[m1(i)], 2);
        dS -= pow(p[p1(i)] - p[i],2);
        dS -= pow(p[i] - p[m1(i)],2);
        dS /= 2*eta;

        double r = exp(-dS);

        if (r > 1) {
            p[i] = yp;
            ++accepted;
        }
        else {
            x = ran.doub();
            if (x < r) {
                p[i] = yp;
                ++accepted;
            }
        }
    }
    return accepted;
}

int PathSim::TailorStep(double eta) {
    unsigned int start = 0;
    unsigned int stop = p1(start);
    bool found = false;
    while (!found && stop < (size - 1))  {
        if (std::abs(p[stop] - p[start] - 0.5) < 0.2*eta) {
            found = true;
        } else {
            stop = p1(stop);
        }
    }

    if (found) {
        // Metropolis test
        S1 stop_p (2*p[start] - p[stop]);
        double dS = pow(p[p1(stop)] - stop_p, 2);
        dS -= pow(p[p1(stop)] - p[stop],2);
        dS /= 2*eta;

        double r = exp(-dS);

        if (r > 1) {
            unsigned int i = start;
            while (i != stop) {
                i = p1(i);
                p[i] = 2*p[start] - p[i];
            }
            return 1;
        }
        else {
            double x = ran.doub();
            if (x < r) {
                unsigned int i = start;
                while (i != stop) {
                    i = p1(i);
                    p[i] = 2*p[start] - p[i];
                }
                return 1;
            }
        }
    }
    return 0;
}
