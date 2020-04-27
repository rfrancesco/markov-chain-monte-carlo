CPP=g++
CPPARGS=-O2 -march=native -std=c++17 -Wpedantic

ran2_class.o : ../../rng/ran2_class.cpp
	$(CPP) $(CPPARGS) -c ../../rng/ran2_class.cpp

s1.o: s1.cpp
	$(CPP) $(CPPARGS) -c s1.cpp 

pathsim.o: pathsim.cpp
	$(CPP) $(CPPARGS) -c pathsim.cpp

local.o: local.cpp
	$(CPP) $(CPPARGS) -c local.cpp

local: local.o s1.o pathsim.o ran2_class.o
	$(CPP) -o xlocal local.o s1.o pathsim.o ran2_class.o

tailor.o : tailor.cpp
	$(CPP) $(CPPARGS) -c tailor.cpp

tailor: tailor.o s1.o pathsim.o ran2_class.o
	$(CPP) -o xtailor tailor.o s1.o pathsim.o ran2_class.o

slab.o: slab.cpp
	$(CPP) $(CPPARGS) -c slab.cpp

slab: slab.o s1.o pathsim.o ran2_class.o
	$(CPP) -o xslab slab.o s1.o pathsim.o ran2_class.o

test: s1.o pathsim.o ran2_class.o
	$(CPP) $(CPPARGS) -c test.cpp
	$(CPP) -o xtest test.o s1.o pathsim.o ran2_class.o
