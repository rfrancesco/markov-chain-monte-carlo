CPP=g++
CPPARGS=-O2 -march=native -std=c++17 -Wpedantic

s1.o: s1.cpp
	$(CPP) $(CPPARGS) -c s1.cpp 

local.o: local.cpp
	$(CPP) $(CPPARGS) -c local.cpp

circle: circle.o s1.o
	$(CPP) -o xcircle local.o s1.o

tailor.o : tailor.cpp
	$(CPP) $(CPPARGS) -c tailor.cpp

tailor: tailor.o s1.o
	$(CPP) -o xtailor tailor.o s1.o
