# QM simulations based on Metropolis algorithm

- harmonic.c: Quantum harmonic oscillator.
  - TODO
    - [ ] Prev/Next lookup arrays
    - [x] Measure correlation functions <x(0)x(τ)>
- circle.cpp: Free particle on a 1D circle with simple local Metropolis algorithm. 
  - This algorithm is not ideal, as it presents an exponential critical slowing down for η -> 0 (issues with ergodicity)
  - Still work in progress
  - TODO: Implement non-local algorithms, maybe parallel tempering.
