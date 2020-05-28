CPP=g++
CPPARGS=-O2 -march=native -std=c++17 -Wpedantic

ran2_class.o : ../rng/ran2_class.cpp
	$(CPP) $(CPPARGS) -c ../rng/ran2_class.cpp

scalar1d.o: scalar1d.cpp
	$(CPP) $(CPPARGS) -c scalar1d.cpp 

thermo.o : thermo.cpp
	$(CPP) $(CPPARGS) -c thermo.cpp

spectrum.o : spectrum.cpp
	$(CPP) $(CPPARGS) -c spectrum.cpp
	
xthermo: thermo.o scalar1d.o ran2_class.o
	$(CPP) -o xthermo thermo.o scalar1d.o ran2_class.o

xspectrum: spectrum.o scalar1d.o ran2_class.o
	$(CPP) -o xspectrum spectrum.o scalar1d.o ran2_class.o

