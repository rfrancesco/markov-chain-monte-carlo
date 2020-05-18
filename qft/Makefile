CPP=g++
CPPARGS=-O2 -march=native -std=c++17 -Wpedantic

ran2_class.o : ../rng/ran2_class.cpp
	$(CPP) $(CPPARGS) -c ../rng/ran2_class.cpp

scalar1d.o: scalar1d.cpp
	$(CPP) $(CPPARGS) -c scalar1d.cpp 

main.o : main.cpp
	$(CPP) $(CPPARGS) -c main.cpp
	
xscalar: main.o scalar1d.o ran2_class.o
	$(CPP) -o xscalar main.o scalar1d.o ran2_class.o

