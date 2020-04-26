CPP=g++
CPPARGS=-O2 -march=native -std=c++17 -Wpedantic

s1.o: s1.cpp
	$(CPP) $(CPPARGS) -c s1.cpp 

pathsim.o: pathsim.o
	$(CPP) $(CPPARGS) -c pathsim.cpp

local.o: local.cpp
	$(CPP) $(CPPARGS) -c local.cpp

local: local.o s1.o
	$(CPP) -o xlocal local.o s1.o pathsim.o

tailor.o : tailor.cpp
	$(CPP) $(CPPARGS) -c tailor.cpp

tailor: tailor.o s1.o
	$(CPP) -o xtailor tailor.o s1.o pathsim.o

test: s1.o pathsim.o
	$(CPP) $(CPPARGS) -c test.cpp
	$(CPP) -o xtest test.o s1.o pathsim.o
