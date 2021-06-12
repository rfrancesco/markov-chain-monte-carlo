CPPC=g++
CPPFLAGS=-O2 -march=native
CC=gcc
CFLAGS=-O2 -march=native
LFLAGS= -lm

xppising: ising.cpp
	$(CPPC) $(CPPFLAGS) ising.cpp -o xppising

xppis2: is2.cpp
	$(CPPC) $(CPPFLAGS) is2.cpp -o xppis2

xising: ising.c
	$(CC) $(CFLAGS) ising.c -o xising $(LFLAGS)

xpar_ising: ising_openmp.cpp
	$(CPPC) $(CPPFLAGS) ising_openmp.cpp -o xpar_ising

potts: potts.cpp
	$(CPPC) $(CPPFLAGS) potts.cpp -o potts

2-potts: 2-potts.cpp
	$(CPPC) $(CPPFLAGS) 2-potts.cpp -o 2-potts
