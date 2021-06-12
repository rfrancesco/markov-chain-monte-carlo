# Classical Ising simulation based on Metropolis algorithm

- ising.c: C implementation
- ising.cpp: C++ implementation (without doubt, much messier than ising.c)
- run.sh: Simple launcher script with GNU Parallel

- potts.cpp: C++ implementation of a N-valued Potts model
- 2-potts.cpp: C++ implementation of two N-valued Potts model, with a small coupling dH = - g * \sum delta(s1[i, j], s2[i, j])
