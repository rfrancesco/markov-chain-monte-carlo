CC=gcc
CCARGS=-O2 -march=native
LARGS=-lm
CPP=g++
CPPARGS=-O2 -march=native -std=c++17 -Wpedantic

harmonic: harmonic.c
	$(CC) $(CCARGS) harmonic.c $(LARGS) -o xharmonic

s1.o: s1.cpp
	$(CPP) $(CPPARGS) -c s1.cpp 

circle.o: circle.cpp
	$(CPP) $(CPPARGS) -c circle.cpp

circle: circle.o s1.o
	$(CPP) -o xcircle circle.o s1.o

circle2.o : circle2.cpp
	$(CPP) $(CPPARGS) -c circle2.cpp

circle2: circle2.o s1.o
	$(CPP) -o xtailor circle2.o s1.o
