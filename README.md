# Markov Chain Monte Carlo
Implementation of Markov Chain Monte Carlo algorithms for Statistical Physics,
Quantum Mechanics, QFT.

The Python analysis scripts included are unorganized and incomplete. 
I will eventually clean them up.

Reports (in Italian) on data generated with this code, and analyzed with Python, are included under `*/relazione/relazione_*.(tex|pdf)`.

Plots are included in `*/relazione/figures/`.

## Contents
- analysis/  
  (Python) Library for data analysis of autocorrelated data
    - analysis.py: Blocking and Bootstrap algorithm with binning with Numpy 
    - linefit_library.py: Simple analytical linear fit (algorithm from A. G. Frodesen, O. Skjeggestad, H. Tøfte, Probability and Statistics in Particle Physics)
    - blocking_complex.py: Blocking analysis for |O_i|, with O_i complex values.
- gaussian/  
  (C) Gaussian number generator based on the Metropolis algorithm
- ising/   
  (C/C++) 2D Ising model simulation (Metropolis algorithm)
  (C++) 2D N-valued Potts model, and coupled Potts lattices (Zn x Zn).
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
