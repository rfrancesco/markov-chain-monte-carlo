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

    const double Neta = atof(argv[1]);   // è essenzialmente la temperatura: Nη = β/mL^2
    const unsigned int N = atoi(argv[2]);
    const double eta = Neta/N;
    const unsigned long int n_measures = atoi(argv[3]);
    const unsigned long int n_skip = atoi(argv[4]);

    const double delta = sqrt(eta);// 0.5; //sqrt(eta);

    unsigned long int m_accepted = 0; // local metropolis acceptance counter
    unsigned long int t_accepted = 0; // nonlocal tailor move acceptance counter

    cout << "#Simulation with Νη = β/mL^2 = " << Neta << ", η = " << eta << ", N = " << N << endl;
    cout << "#Taking " << n_measures << " measures every " << n_skip << "sweeps + nonlocal move" << endl;
    cout << "#Local Metropolis algorithm (delta = " << delta << " ) + nonlocal tailor move each sweep" << endl;
    cout << "#C1(τ) = <cos(2πτ)cos(0)>_c, C2(τ) = <cos(4πτ)cos(0)>_c" << endl;
    cout << "#Q\t[C1(0...)]\t[C2(0...)]" << endl;

    PathSim* p = new PathSim(N);

    for (unsigned int i = 0; i < n_measures; ++i) {
        // Markov steps
        for (unsigned int j = 0; j < n_skip; ++j) {
            m_accepted += p->MetropolisSweep(eta, delta);
            t_accepted += p->TailorStep(eta);
        }
        // Printing observables
        //printf("%.3f\t", p->winding()); // winding number
        for (unsigned int k = 0; k < N/(Neta); ++k)
            cout << p->cos_correlator(k) << "\t";
        double cos2_avg2 = pow(p->cos2_avg(),2);
        for (unsigned int k = 0; k < N/(Neta); ++k)
            cout << p->cos2_correlator(k) - cos2_avg2 << "\t";
        cout << "\n";
    }



    const unsigned long int m_steps = n_measures*N*n_skip;
    const unsigned long int t_steps = n_measures*n_skip;
    const unsigned long int tot_steps = m_steps + t_steps;
    const unsigned long int tot_accepted = m_accepted + t_accepted;


    cout << "#Local Metropolis acceptance: " << m_accepted << " out of " << m_steps << " (" << static_cast<double>(m_accepted)/m_steps << ")" << endl;
    cout << "#Nonlocal move acceptance: " << t_accepted << " out of " << t_steps << " (" << static_cast<double>(t_accepted)/t_steps << ")" << endl;
    cout << "#Total acceptance: " << tot_accepted << " out of " << tot_steps << " (" << static_cast<double>(tot_accepted)/tot_steps << ")" << endl;

    return 0;
}
