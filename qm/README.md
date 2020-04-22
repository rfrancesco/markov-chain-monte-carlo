# QM simulations based on Metropolis algorithm


- harmonic.c: Quantum harmonic oscillator.
  - TODO
    - [ ] Prev/Next lookup arrays
    - [x] Measure correlation functions <x(0)x(τ)>
    
    
- circle.cpp: Free particle on a 1D circle with simple local Metropolis algorithm. 
  This algorithm is not ideal, as it presents an exponential critical slowing down for η -> 0 (issues with ergodicity)
- circle2.cpp: Free particle on a 1D circle with local Metropolis + nonlocal "Tailor" move
  Reference: C. Bonati, M. D'Elia, _Topological critical slowing down:
    Variations on a toy model_, Phys. Rev. E 98, 013308 (2018)

- TODO:
  - [ ] Implement Parallel Tempering algorithm
    - [ ] ^- prerequisite: write cleaner code for the Path object. It does work,
        but the structure of the program could be more elegant.
    - [ ] ^- probably an interesting exercise to learn the basics of multiprocessing
    
  
- s1.hpp/cpp: C++ data type implementing the circle group S1
  - x in [-1/2, 1/2)
  - Addition respects group composition
  - Group inverse implemented by unary minus
  - Subtraction naturally implements distance with sign
  - Casts naturally to a double when needed
  - Still work in progress (probably needs const correctness, etc.)
  
