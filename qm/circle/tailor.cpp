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

/*  Alcune note sulla classe S1:
 *  ===============================================================
 *  S1 x(0.3) crea un oggetto S1 con valore 0.3.
 *  S1 rappresenta gli elementi del gruppo S_1 come double compresi
 *  tra -1/2 e 1/2.
 *
 *  Costruzione: explicit S1 a(double);
 *  !!! Possiamo scrivere S1 a; a = 0.2; ma non S1 a = 0.2; !!!
 *  oppure S1 a; <- viene automaticamente assegnato a 0
 *
 *  Operazioni:
 *  + : composizione di gruppo (+(S1, S1), +(S1, double), +(double, S1))
 *  - : distanza con segno (-(S1, double))
 *  - : inverso nel gruppo (-(S1)), coerente con la sottrazione
 *  = : assegnazione (anche cast double -> S1)
 *  Il risultato è un altro S1.
 *
 *  stream << : es. cout << a; comportamento consistente con gli altri
 *              tipi numerici, senza cast implicito a double.
 *
 *  Può essere castato a double, implicitamente o esplicitamente con
 *  static_cast<double>(a).
 *
 *  Warning: il risultato delle operazioni è di tipo S1, quindi ci vuole
 *  cautela nel maneggiare S1 e double: ad esempio:
 *
 *  Ok: ==========================
 *  S1* p = new S1[n];
 *  // calcolo del numero di avvolgimenti
 *  double d = 0;
 *  for (i...)
 *      d += p[i+1] - p[i] <- il risultato della sottrazione viene castato a double
 *                            la somma += effettuata è quella di double
 *
 *  cout << d; <- Esempio di risultato: "5"
 *
 *  No: ==========================
 *  S1* p = new S1[3];
 *
 *  p[0] = 0;
 *  p[1] = 0.4;
 *  p[2] = -0.4;
 *
 *  // calcolo del numero di avvolgimenti
 *  double d = (p[2] - p[1]) + (p[1] - p[0]) + (p[0] - p[2])
 *                           ^- somma S1 +(S1, S1)
 *  atteso?   d = (double) 0.2 + (double) 0.4 + (double) 0.4 = 1
 *  ottenuto? d = (S1) 0.2 + (S1) 0.4 + (S1) 0.4 = (S1) 1 ~ (S1) 0 -> (double) 0
 *
 *  Questo perché gli elementi nella parentesi sono S1 e sommano come S1,
 *  non come double. Per ovviare si può o sommarli separatamente (come nel codice OK)
 *  oppure castare ogni parentesi a double, che però è pesante e poco leggibile.
 *
 *
 *  double d = static_cast<double>(p[2] - p[1]) + static_cast<double>(p[1] - p[0]) + static_cast<double>(p[0] - p[2])
 *
 *  Un altro (e ultimo) caveat: in generale, (x - y)^2 != x^2 + y^2 - 2xy
 *  La moltiplicazione non ha un significato intrinseco in S1, e la distanza non gode di questa proprietà.
 *  */

// Let us initialize our RNG with a random seed
// Using <ctime> time(NULL); would not be enough, as simulations
// that are run in parallel with GNU Parallel could have the same seed

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
        printf("%.3f\t", p->winding()); // winding number
        /*for (unsigned int k = 0; k < 1.5*N/(10*Neta); ++k)
            cout << p->cos_correlator(k) << "\t";
        double cos2_avg2 = pow(p->cos2_avg(),2);
        for (unsigned int k = 0; k < 1.5*N/(10*Neta); ++k)
            cout << p->cos2_correlator(k) - cos2_avg2 << "\t";*/
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
