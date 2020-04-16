CPPC=g++
CPPFLAGS=-O2 -march=native
CC=gcc
CFLAGS=-O2 -march=native
LFLAGS= -lm

xppising: ising.cpp
	$(CPPC) $(CPPFLAGS) ising.cpp -o xppising

xising: ising.c
	$(CC) $(CFLAGS) ising.c -o xising $(LFLAGS)
