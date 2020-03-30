CPPC=g++
CPPFLAGS=-O2 -march=native
CC=gcc
CFLAGS=-O2 -march=native -lm

ising: ising.cpp
	$(CPPC) $(CPPFLAGS) ising.cpp -o xising

cising: ising.c
	$(CC) $(CFLAGS) ising.c -o xcising
