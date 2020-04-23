#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <cmath>

#include "s1.hpp"
#include "../../rng/ran2_double.cpp"

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
        // Observables
        double winding() const;
        double cos_avg() const;
        double cos_correlator(unsigned int) const;
        double cos2_avg() const;
        double cos2_correlator(unsigned int) const;
        // Markov steps
        int MetropolisSweep(double, double);
        int TailorStep(double);
};

double Path::cos_avg() const {
    double result = 0;
    // Calculate <x^2> as average over entire path
    for (unsigned int i = 0; i < size; ++i)
        result += cos(2*M_PI*p[i]);
    result /= size;
    return result;
}

double Path::cos_correlator(unsigned int k) const {
    double result = 0;
    // Calculate C(k) as average of C(0, k), C(1, 1+k) ... C(s-1, s-1+k)
    for (unsigned int i = 0; i < size; ++i) {
        unsigned int ik = i;
        for (unsigned int j = 0; j < k; ++j)
            ik = p1(ik);
        result += cos(2*M_PI*p[ik])*cos(2*M_PI*p[i]);
    }
    result /= size;
    return result;
}

double Path::cos2_avg() const {
    double result = 0;
    // Calculate <x^2> as average over entire path
    for (unsigned int i = 0; i < size; ++i)
        result += cos(4*M_PI*p[i]);
    result /= size;
    return result;
}

double Path::cos2_correlator(unsigned int k) const {
    double result = 0;
    // Calculate C(k) as average of C(0, k), C(1, 1+k) ... C(s-1, s-1+k)
    for (unsigned int i = 0; i < size; ++i) {
        unsigned int ik = i;
        for (unsigned int j = 0; j < k; ++j)
            ik = p1(ik);
        result += cos(4*M_PI*p[ik])*cos(4*M_PI*p[i]);
    }
    result /= size;
    return result;
}

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
    printf("%.3f\t", winding());
    double cos_avg2 = pow(cos_avg(),2);
    for (unsigned int i = 0; i < size; i++)
        cout << cos_correlator(i) - cos_avg2 << "\t"; // connected correlator
    double cos2_avg2 = pow(cos2_avg(),2);
    for (unsigned int i = 0; i < size; i++)
        cout << cos2_correlator(i) - cos2_avg2 << "\t"; // connected correlator
    cout << '\n';
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
        dS /= 2*eta;

        double r = exp(-dS);

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

int Path::TailorStep(double eta) {
    //cout << "Attempting tailor (" << winding() << ")...";
    //unsigned int start = floor(size * ran2());
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
            double x = ran2();
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
    //cout << "failed" << endl;
    return 0;
}

int main(int argc, char** argv) {

    if (argc != 5) {
        cout << "Usage: 4 arguments: Neta, N, n_measures, n_skip" << endl;
        return EXIT_FAILURE;
    }

    ran2_init();
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
    cout << "#Q\t";
    for (unsigned int i = 0; i < N; ++i)
        cout << "C1(" << i << ")\t";
    for (unsigned int i = 0; i < N; ++i)
        cout << "C2(" << i << ")\t";
    cout << endl;

    Path* p = new Path(N);

    for (unsigned int i = 0; i < n_measures; ++i) {
        for (unsigned int j = 0; j < n_skip; ++j) {
            m_accepted += p->MetropolisSweep(eta, delta);
            t_accepted += p->TailorStep(eta);
        }
        p->print();
    }

    const unsigned long int m_steps = n_measures*N*n_skip;
    const unsigned long int t_steps = n_measures*n_skip;
    const unsigned long int tot_steps = m_steps + t_steps;
    const unsigned long int tot_accepted = m_accepted + t_accepted;


    cout << "#Local Metropolis acceptance: " << m_accepted << " out of " << m_steps << " (" << static_cast<double>(m_accepted)/m_steps << ")" << endl;
    cout << "#Nonlocal move acceptance: " << t_accepted << " out of " << t_steps << " (" << static_cast<double>(t_accepted)/t_steps << ")" << endl;
    cout << "#Total acceptance: " << tot_accepted << " out of " << tot_steps << " (" << static_cast<double>(tot_accepted)/tot_steps << ")" << endl;

    ran2_save();
    return 0;
}
