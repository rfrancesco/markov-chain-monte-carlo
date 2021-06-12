#include <iostream>
#include <cmath>
#include <complex>
#include "../rng/ran2_double.cpp"

const unsigned int N_POTTS = 3;


/* int random_spin() { 		// Genera casualmente \pm 1 con ran2
	double x = ran2();
	if (x > 0.5)
		return 1;
	else
		return -1;
} */

unsigned int random_spin() {
	unsigned int x = floor(N_POTTS*ran2());
	return (x % N_POTTS);
}

int delta(int x, int y) {
	if (x == y) {
		return 1;
	} else {
		return 0;
	}
}


class Lattice {
	private:
		unsigned int **matrix;			// Matrice di interi
		unsigned int size;		// Dimensioni della matrice
		unsigned int **geometry;	// Implementa le condizioni periodiche
	public:
		Lattice (unsigned int, int);		// Costruttore
		~Lattice();			// Distruttore
		const unsigned int& operator() (int, int) const;	// Operatore lattice(i,j) (lettura)
		unsigned int& operator() (int, int);	// Operatore lattice(i,j) (scrittura)
		unsigned int p1(int) const;		// operatore +1 size-periodico
		unsigned int m1(int) const;		// operatore -1 size-periodico
};


const unsigned int& Lattice::operator()(int i, int j) const {
	return matrix[i][j];
}

unsigned int& Lattice::operator()(int i, int j) {
	return matrix[i][j];
}

unsigned int Lattice::p1(int i) const {
	return geometry[0][i];
}

unsigned int Lattice::m1(int i) const {
	return geometry[1][i];
}

Lattice::Lattice(unsigned int n, int init) {
	// Constructor dell'oggetto Lattice
	// Crea un cammino di dimensione s e lo inizializza a 0
	size = n;
	matrix = nullptr;
	geometry = nullptr;
	if (n != 0) {
		// Alloco la matrice
		matrix = new unsigned int*[n];
		for (unsigned int i = 0; i < n; i++)
			matrix[i] = new unsigned int[n];
		// Alloco la matrice interna delle condizioni al bordo
		geometry = new unsigned int*[2];
		geometry[0] = new unsigned int[n];
		geometry[1] = new unsigned int[n];
		for (unsigned int i = 0; i < n; i++) {
			geometry[0][i] = i + 1;
			geometry[1][i] = i - 1;
		}
		geometry[0][n-1] = 0;
		geometry[1][0] = n - 1;
	}
	// Inizializzazione: 0 a freddo, 1 a caldo, altrimenti errore.
	// Andrebbe spostato in Simulation, ma funziona uguale
	// Magari ci lascio un parametro qui in futuro, per riciclare il codice
	switch (init) {
		case 0:		for (unsigned int i = 0; i < n; i++)
					for (unsigned int j = 0; j < n; j++)
						matrix[i][j] = 0;
				std::cout << "#Lattice initialized with 0\n";
				break;
		case 1: 	for (unsigned int i = 0; i < n; i++)
					for (unsigned int j = 0; j < n; j++)
						matrix[i][j] = random_spin();
				std::cout << "#Lattice initialized with random spins\n";
				break;
		default: 	std::cerr << "Illegal value in lattice init. Exiting..." << std::endl;
				exit(EXIT_FAILURE);
	}
}


Lattice::~Lattice() {
	// Distruttore dell'oggetto Lattice
	// Dealloca la matrice
	for (unsigned int i = 0; i < size; i++)
		delete [] matrix[i];
	delete [] matrix;
	// Dealloco la matrice interna delle condizioni al bordo
	delete [] geometry[0];
	delete [] geometry[1];
	delete [] geometry;
}

class Simulation { 			// L'oggetto che gestisce la simulazione
	private:
		unsigned int size;		// grandezza della matrice
		unsigned int volume;		// volume della matrice
		unsigned long long accepted;		// passi accettati
		unsigned long long steps;		// passi fatti
		double beta;		// parametro beta
		double extfield;	// Campo magnetico esterno
		double g;		// l1-l2 coupling
		Lattice *l1, *l2;		// Lattice
		// Math auxiliary stuff
		double delta[N_POTTS][N_POTTS]; // delta function
		std::complex<double> roots[N_POTTS];
	public:
		Simulation(unsigned int, double, double, double, int);	// Costruttore
		~Simulation();				// Distruttore
		void Metrostep_lattice(unsigned int lattice);
		void Metrostep();			// Singolo passo Metropolis
		double Acceptance() const;			// Restituisce l'accettanza
		double Energy() const ;			// Misura l'energia
		void Magnetization() const;			// Stampa le magnetizzazioni
		void MixedMagnetization() const;
};

Simulation::Simulation(unsigned int n, double b, double h, double g_c, int init) {
	// Costruttore dell'oggetto Simulation
	// Importa tutti i parametri necessari, e costruisce un oggetto Lattice
	size = n;
	extfield = h;
	volume = size*size;
	beta = b;
	g = g_c;
	steps = 0;
	accepted = 0;
	l1 = new Lattice(size, init);
	l2 = new Lattice(size, init);

	// Populate delta[][] function
	for (unsigned int a = 0; a < N_POTTS; a++) {
		for (unsigned int b = 0; b < N_POTTS; b++) {
			if (a == b)
				delta[a][b] = 1;
			else
				delta[a][b] = 0;
		}
	}

	// Populate roots[]
	const double k = 2*M_PI/N_POTTS;
	const std::complex<double> i(0,1);
	for (unsigned int a = 0; a < N_POTTS; a++)
		roots[a] = exp(i*(k*a));
}

Simulation::~Simulation() {
	// Distruttore dell'oggetto Simulation
	// Dealloca l'oggetto Lattice
	delete l1;
	delete l2;
}

double Simulation::Acceptance() const {
	double a = static_cast<double>(accepted)/steps;
	return a;
}

void Simulation::Metrostep_lattice(unsigned int lattice) {
	Lattice *l, *other;
	switch (lattice) {
		case 0: l = l1;
			other = l2;
			break;
		case 1: l = l2;
			other = l1;
			break;
		default: exit(EXIT_FAILURE);
	}

	// Scelgo un sito casualmente
	int i = floor((ran2() * size));
	int j = floor((ran2() * size));

	// Calcolo ++i, --i, ++j, --j
	const int jp = l->p1(j);
	const int jm = l->m1(j);
	const int ip = l->p1(i);
	const int im = l->m1(i); 
	const unsigned int old_value = (*l)(i,j);
	const unsigned int new_value = random_spin();
	const unsigned int other_value = (*other)(i,j);

	// Calcolo i contributi all'energia
	double dH = - (delta[new_value][(*l)(i, jp)] + delta[new_value][(*l)(i, jm)] + delta[new_value][(*l)(ip, j)] + delta[new_value][(*l)(im, j)]) 
	            + (delta[old_value][(*l)(i, jp)] + delta[old_value][(*l)(i, jm)] + delta[old_value][(*l)(ip, j)] + delta[old_value][(*l)(im, j)]);
	double dH_extfield = - extfield * (delta[0][new_value] - delta[0][old_value]);
	double dH_g = - g * (delta[other_value][new_value] - delta[other_value][new_value]);

	double dbH = beta * (dH + dH_extfield + dH_g);

	if (dbH < 0) {
		(*l)(i,j) = new_value;
		accepted++;
	} else {
		double r = exp(-dbH);
		double x = ran2();
		if (x < r) {
			(*l)(i,j) = new_value;
			accepted++;	// aggiorno l'accettanza
		}
	}
	steps++;		// aggiorno il numero di passi
}

void Simulation::Metrostep() {
	Metrostep_lattice(0);
	Metrostep_lattice(1);
}

double Simulation::Energy() const {
	double e = 0;
	double f;
	for (unsigned int i = 0; i < size; i++) {
		for (unsigned int j = 0; j < size; j++) {
			// lattice 1
			int jp = l1->p1(j);
			int jm = l1->m1(j);
			int ip = l1->p1(i);
			int im = l1->m1(i); 
			unsigned int s1 = (*l1)(i,j);
			double dE = - (delta[s1][(*l1)(i, jp)] + delta[s1][(*l1)(i, jm)] + delta[s1][(*l1)(ip, j)] + delta[s1][(*l1)(im, j)])
			     - extfield*delta[0][s1];
			// lattice 2
			jp = l2->p1(j);
			jm = l2->m1(j);
			ip = l2->p1(i);
			im = l2->m1(i); 
			unsigned int s2 = (*l2)(i,j);
			dE += - (delta[s2][(*l2)(i, jp)] + delta[s2][(*l2)(i, jm)] + delta[s2][(*l2)(ip, j)] + delta[s2][(*l2)(im, j)])
			     - extfield*delta[0][s2];
			dE /= 2;
			// interaction
			dE += - beta*g*delta[s1][s2];
			e += dE;	
		}
	}
	e /= volume;
	return e;
}

void Simulation::Magnetization() const {
	double m1[N_POTTS], m2[N_POTTS];
	for (unsigned int a = 0; a < N_POTTS; a++) {
		m1[a] = 0;
		m2[a] = 0;
	}
	for (unsigned int i = 0; i < size; i++) {
		for (unsigned int j = 0; j < size; j++) {
			m1[(*l1)(i,j)]++;				
			m2[(*l2)(i,j)]++;				
		}
	}
	std::complex<double> mag1 = 0;
	std::complex<double> mag2 = 0;
	for (unsigned int a = 0; a < N_POTTS; a++) {
		m1[a] /= volume;
		m2[a] /= volume;
		mag1 += m1[a]*roots[a];
		mag2 += m2[a]*roots[a];
	}
	std::cout << std::real(mag1) << "\t" << std::imag(mag1) << "\t";
	std::cout << std::real(mag2) << "\t" << std::imag(mag2) << "\t";
}

unsigned int mod(int i, unsigned int n) {
	while (i < 0)
		i += n;
	while (i > n) 
		i -= n;
	return i;
}

void Simulation::MixedMagnetization() const {
	double mm[N_POTTS], mmd[N_POTTS];
	for (unsigned int a = 0; a < N_POTTS; a++) {
		mm[a] = 0;
		mmd[a] = 0;
	}
	for (unsigned int i = 0; i < size; i++) {
		for (unsigned int j = 0; j < size; j++) {
			unsigned int mm_ = mod( (*l1)(i,j) + (*l2)(i,j), N_POTTS);
			unsigned int mmd_ = mod( (*l1)(i,j) - (*l2)(i,j), N_POTTS);
			mm[mm_]++;				
			mmd[mmd_]++;				
		}
	}
	std::complex<double> mag1 = 0;
	std::complex<double> mag2 = 0;
	for (unsigned int a = 0; a < N_POTTS; a++) {
		mm[a] /= volume;
		mmd[a] /= volume;
		mag1 += mm[a]*roots[a];
		mag2 += mmd[a]*roots[a];
	}
	std::cout << std::real(mag1) << "\t" << std::imag(mag1) << "\t";
	std::cout << std::real(mag2) << "\t" << std::imag(mag2);
}


int main(int argc, char** argv) {

	if (argc != 8) {
		std::cout << "Usage: 7 arguments: n_measures, n_skip, n, beta, extfield, g, init" << std::endl;
		std::cout << "n_skip * N^2: metropolis steps between two measures" << std::endl;
		std::cout << "init = (0: cold, 1: hot)" << std::endl;
		return EXIT_FAILURE;
	}

	// Inizializzo ran2
	ran2_init();
	
	// Stampa i double con la stessa precisione
	// con cui vengono immagazzinati in memoria
	std::cout.precision(17); 
	
	// Parametri. Da prendere come argomenti
	// in modo da non dover ricompilare
	const unsigned int n_measures = atoi(argv[1]);
	const unsigned int n = atoi(argv[3]);
	const unsigned int n_skip = atoi(argv[2]) * n * n;
	const double beta = atof(argv[4]);
	const double extfield = atof(argv[5]);
	const double g = atof(argv[6]);
	const int init = atoi(argv[7]);	// 0: freddo, 1: caldo
	// Inizializzo la simulazione
	Simulation sim(n, beta, extfield, g, init);

	// Inizio misure
	std::cout << "#Beginning simulation with beta = " << beta << ", N = " << n << ", extfield = " << extfield << ", g = " << g << " \n";
	std::cout << "#Calculating " << n_measures << "measures, skipping " << n_skip << "cycles per measure." << std::endl;
	std::cout << "#Energy\tMagnetization\n";

	for (unsigned int i = 0; i < n_measures; i++) {
		for (unsigned int j = 0; j < n_skip; j++)
			sim.Metrostep();
		std::cout << sim.Energy() << "\t";
		sim.Magnetization();
		sim.MixedMagnetization();
		std::cout << std::endl;
		std::cout << std::flush;
	} 

	// Accettanza
	std::cout << "#Generated " << n_measures << " measures with acceptance: " << sim.Acceptance() << std::endl;

	// Salvo stato del generatore
	ran2_save(); 
	return EXIT_SUCCESS;
}

