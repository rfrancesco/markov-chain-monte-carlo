# Markov Chain Monte Carlo
Implementation of Markov Chain Monte Carlo algorithms for Statistical Physics,
Quantum Mechanics, QFT (work in progress)

## Contents
- analysis/  
  (Python) Library for data analysis of autocorrelated data
    - analysis.py: Blocking and Bootstrap algorithm with binning with Numpy 
    - linefit_library.py: Simple analytical linear fit (algorithm from A. G. Frodesen, O.
   Skjeggestad, H. Tøfte, Probability and Statistics in Particle Physics)
- gaussian/  
  (C) Gaussian number generator based on the Metropolis algorithm
  - gauss.c: C implementation
- ising/   
  (C/C++) 2D Ising model simulation (Metropolis algorithm)
    - ising.c: C implementation
    - ising.cpp: C++ implementation
    - run.sh: Very simple launcher with GNU Parallel
- qm/  
  (C) QM simulations
    - harmonic.c: Quantum harmonic oscillator
    - s1.hpp/cpp: C++ data type implementing the circle group S1
    - circle.cpp: Quantum free particle on a circle (S1), local Metropolis
    - circle2.cpp: Quantum free particle on a circle (S1), local Metropolis +
      nonlocal Tailor move
        
The random number generator ran2() (rng/ran2.c, rng/ran2_double.cpp) is not included due to licensing concerns.
