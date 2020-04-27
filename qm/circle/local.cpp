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
// This is not good if you want to _know_ the state of the RNG
// but it is certainly better than saving/loading from a file,
// which could ruin parallel simulations

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
