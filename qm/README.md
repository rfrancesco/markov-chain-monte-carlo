# QM simulations based on Metropolis algorithm

- harmonic.c: Quantum harmonic oscillator.
  - TODO
    - [ ] Prev/Next lookup arrays
    - [x] Measure correlation functions <x(0)x(τ)>
- circle.cpp: Free particle on a 1D circle with simple local Metropolis algorithm. 
  - This algorithm is not ideal, as it presents an exponential critical slowing down for η -> 0 (issues with ergodicity)
  - Still work in progress
  - TODO: Implement non-local algorithms, maybe parallel tempering.
- s1.hpp/cpp: C++ data type implementing the circle group S1
  - x in [-1/2, 1/2)
  - Addition respects group composition
  - Inverse
  - Subtraction implements distance with sign
  - Casts naturally to a double when needed
  - Still work in progress (probably needs const correctness, etc.)
  
