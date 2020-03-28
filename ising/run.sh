#!/bin/bash

# Run simulations with GNU Parallel

BINARY=xising
N_JOBS=2	# How many jobs should be run?
NS="10 20 30 40 50"	# lattice size
BETAS="0.3 0.4 0.5"	# values of beta to simulate


parallel -j $N_JOBS --progress --results "{1}/ising_{2}" ./$BINARY 100000 {1} 20 {2} 0 1 ::: $NS ::: $BETAS
