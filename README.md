# Markov Chain Monte Carlo
Implementation of Markov Chain Monte Carlo algorithms for Statistical Physics, Quantum Mechanics, QFT

## Contents
- analysis/: (Python) Library for data analysis of autocorrelated data
    - analysis.py: Blocking and Bootstrap algorithm with binning with Numpy 
    - linefit_library.py: Simple analytical linear fit (algorithm from A. G. Frodesen, O.
   Skjeggestad, H. Tøfte, Probability and Statistics in Particle Physics)
- ising/: (C/C++) 2D Ising model simulation (Metropolis algorithm)
    - ising.c: C implementation
    - ising.cpp: C++ implementation
    - run.sh: Very simple launcher with GNU Parallel
 
The random number generator ran2() (rng/ran2.c, rng/ran2_double.cpp) is not included due to licensing concerns.
