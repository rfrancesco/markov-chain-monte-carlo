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
  (C/C++) QM simulations
    - harmonic.c: Quantum harmonic oscillator
    - circle/: (C++) Quantum free particle on a circle.  
      Contains various algorithms, to study and address the
      exponential critical slowing down for η -> 0, in the β >> 1 regime,
      and to study the spectrum of the Hamiltonian.
- qft/
  (C++) 1+1D Scalar quantum field
  Contains tools to study the thermodynamics and the spectrum of the field,
  by using wall-wall correlators.
        
The random number generator ran2() (rng/ran2.c, rng/ran2_double.cpp) is not included due to licensing concerns.
