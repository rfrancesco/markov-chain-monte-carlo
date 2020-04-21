#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <cmath>

#include "s1.hpp"
#include "../rng/ran2_double.cpp"

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





class Path {
    private:
        S1 * p;
        unsigned int size;
        unsigned int **geometry;
    public:
        Path(unsigned int);
        ~Path();
        const S1& operator[] (unsigned int) const;
        S1& operator[] (unsigned int);
        unsigned int p1(unsigned int) const;
        unsigned int m1(unsigned int) const;
        void print() const;
        double winding() const;
        int MetropolisSweep(double, double);
};

double Path::winding() const {
    double q = 0;
    for (unsigned int i = 0; i < size; i++)
        q += p[p1(i)] - p[i];
    return q;
}


Path::Path(unsigned int s) {
    p = new S1[s];
    // Alloco la matrice interna delle condizioni al bordo
    geometry = new unsigned int*[2];
    geometry[0] = new unsigned int[s];
    geometry[1] = new unsigned int[s];
    for (unsigned int i = 0; i < s; i++) {
        geometry[0][i] = i + 1;
        geometry[1][i] = i - 1;
    }
    geometry[0][s-1] = 0;
    geometry[1][0] = s - 1;
    size = s;
    for (unsigned int i = 0; i < size; i++)
        p[i] = ran2();
}

Path::~Path() {
    delete [] geometry[0];
    delete [] geometry[1];
    delete [] p;
}

void Path::print() const {
    printf("%.3f\n", winding());
}


unsigned int Path::p1(unsigned int i) const {
    return geometry[0][i];
}

unsigned int Path::m1(unsigned int i) const {
    return geometry[1][i];
}

const S1& Path::operator[](unsigned int i) const {
    return p[i];
}

S1& Path::operator[](unsigned int i) {
    return p[i];
}

int Path::MetropolisSweep(double eta, double delta) {
    unsigned int accepted = 0;
    for (unsigned int i = 0; i < size; i++) {
        double x = delta*(2*ran2() - 1);
        S1 yp = p[i] + x;

        double dS = pow( p[p1(i)] - yp , 2);
        dS += pow(yp - p[m1(i)], 2);
        dS -= pow(p[p1(i)] - p[i],2);
        dS -= pow(p[i] - p[m1(i)],2);

        //double r = exp(-dS * (2*PI2/eta));
        double r = exp(-dS /(2*eta));

        if (r > 1) {
            p[i] = yp;
            ++accepted;
        }
        else {
            x = ran2();
            if (x < r) {
                p[i] = yp;
                ++accepted;
            }
        }
    }
    return accepted;
}


int main() {
    ran2_init();
    const unsigned int N = 600;
    const double Neta = 5;
    const double eta = Neta/N;   // è essenzialmente la temperatura: Nη = β/mR^2

    const double delta = 0.5;//sqrt(eta);

    const unsigned int n_measures = 2e5;
    unsigned long int accepted = 0;

    cout << "#Simulation with Νη = β/mR^2 = " << Neta << ", η = " << eta << ", N = " << N << endl;

    Path* p = new Path(N);

    for (unsigned int i = 0; i < n_measures; ++i) {
        accepted += p->MetropolisSweep(eta, delta);
        p->print();
    }

    cout << "#Acceptance: " << accepted << "out of " << N*n_measures << " (" << static_cast<double>(accepted)/(N*n_measures) << ")" << endl;

    ran2_save();
    return 0;
}
