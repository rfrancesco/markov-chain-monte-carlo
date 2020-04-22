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
- ising/   
  (C/C++) 2D Ising model simulation (Metropolis algorithm)
- qm/  
  (C) QM simulations
    - harmonic.c: Quantum harmonic oscillator
    - circle/: Quantum free particle on a circle.  
      Contains various algorithms, implemented in C++, to study and address the
      exponential critical slowing down for η -> 0, in the β >> 1 regime.
        
The random number generator ran2() (rng/ran2.c, rng/ran2_double.cpp) is not included due to licensing concerns.
