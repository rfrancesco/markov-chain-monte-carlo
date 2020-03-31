#include <iostream>
#include <cmath>
#include "../rng/ran2_double.cpp"


int random_spin() { 		// Genera casualmente \pm 1 con ran2
	double x = ran2();
	if (x > 0.5)
		return 1;
	else
		return -1;
}

class Lattice {
	private:
		int **matrix;			// Matrice di interi
		unsigned int size;		// Dimensioni della matrice
		unsigned int **geometry;	// Implementa le condizioni periodiche
	public:
		Lattice (unsigned int, int);		// Costruttore
		~Lattice();			// Distruttore
		const int& operator() (int, int) const;	// Operatore lattice(i,j) (lettura)
		int& operator() (int, int);	// Operatore lattice(i,j) (scrittura)
		unsigned int p1(int) const;		// operatore +1 size-periodico
		unsigned int m1(int) const;		// operatore -1 size-periodico
};


const int& Lattice::operator()(int i, int j) const {
	return matrix[i][j];
}

int& Lattice::operator()(int i, int j) {
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
		matrix = new int*[n];
		for (unsigned int i = 0; i < n; i++)
			matrix[i] = new int[n];
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
						matrix[i][j] = 1;
				std::cout << "#Lattice initialized with 1\n";
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
		Lattice *l;		// Lattice
	public:
		Simulation(unsigned int, double, double, int);	// Costruttore
		~Simulation();				// Distruttore
		void Metrostep();			// Singolo passo Metropolis
		double Acceptance() const;			// Restituisce l'accettanza
		double Energy() const ;			// Misura l'energia
		double Magnetization() const;			// Misura la magnetizzazione
};

Simulation::Simulation(unsigned int n, double b, double h, int init) {
	// Costruttore dell'oggetto Simulation
	// Importa tutti i parametri necessari, e costruisce un oggetto Lattice
	size = n;
	extfield = h;
	volume = size*size;
	beta = b;
	steps = 0;
	accepted = 0;
	l = new Lattice(size, init);
}

Simulation::~Simulation() {
	// Distruttore dell'oggetto Simulation
	// Dealloca l'oggetto Lattice
	delete l;
}

double Simulation::Acceptance() const {
	double a = static_cast<double>(accepted)/steps;
	return a;
}

void Simulation::Metrostep() {
	// Scelgo un sito casualmente
	int i = floor((ran2() * size));
	int j = floor((ran2() * size));

	// Calcolo ++i, --i, ++j, --j
	const int jp = l->p1(j);
	const int jm = l->m1(j);
	const int ip = l->p1(i);
	const int im = l->m1(i); 

	// Calcolo i contributi all'energia
	double f = (*l)(i, jp) + (*l)(i, jm) + (*l)(ip, j) + (*l)(im, j);
	f = beta*(f + extfield);

	double r = exp(-2.0 * (*l)(i,j) * f);

	// Test di Metropolis
	double x = ran2();
	if (x < r) {
		(*l)(i,j) = -(*l)(i,j);
		accepted++;	// aggiorno l'accettanza
	}
	steps++;		// aggiorno il numero di passi
}

double Simulation::Energy() const {
	double e = 0;
	double f;
	for (unsigned int i = 0; i < size; i++) {
		for (unsigned int j = 0; j < size; j++) {
			f = (*l)(i, l->p1(j)) + (*l)(i, l->m1(j)) + (*l)(l->p1(i), j) + (*l)(l->m1(i), j);
			e -= 0.5*f*(*l)(i,j);
			e -= extfield*(*l)(i,j);	
		}
	}
	e /= volume;
	return e;
}

double Simulation::Magnetization() const {
	double m = 0;
	for (unsigned int i = 0; i < size; i++)
		for (unsigned int j = 0; j < size; j++)
			m += (*l)(i,j);
	m /= volume;
	return m;
}




int main(int argc, char** argv) {

	if (argc != 7) {
		std::cout << "Usage: 6 arguments: n_measures, n_skip, n, beta, extfield, init" << std::endl;
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
	const int init = atoi(argv[6]);	// 0: freddo, 1: caldo
	// Inizializzo la simulazione
	Simulation sim(n, beta, extfield, init);

	// Inizio misure
	std::cout << "#Beginning simulation with beta = " << beta << ", N = " << n << ", extfield = " << extfield << " \n";
	std::cout << "#Calculating " << n_measures << "measures, skipping " << n_skip << "cycles per measure." << std::endl;
	std::cout << "#Energy\tMagnetization\n";

	for (unsigned int i = 0; i < n_measures; i++) {
		for (unsigned int j = 0; j < n_skip; j++)
			sim.Metrostep();
		std::cout << sim.Energy() << "\t" << sim.Magnetization() << "\n";
	} 

	// Accettanza
	std::cout << "#Generated " << n_measures << " measures with acceptance: " << sim.Acceptance() << std::endl;

	// Salvo stato del generatore
	ran2_save(); 
	return EXIT_SUCCESS;
}

