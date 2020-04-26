# QM simulations based on Metropolis algorithm


- harmonic.c: Quantum harmonic oscillator.
  - TODO
    - [ ] Prev/Next lookup arrays
    - [x] Measure correlation functions <x(0)x(τ)>
    
   
- circle/: Quantum free particle on a 1D circle  
  Contents:
  - pathsim.hpp/cpp: C++ class for simulating a free particle on a 1D circle with MCMC
    - The engine of these simulations
    - Measurements implemented as public functions
    - MCMC algorithms implemented as public functions
  - s1.hpp/cpp: C++ data type implementing the circle group S1  
    Features:
    - Double precision
    - Wraps any double in the domain [-1/2, +1/2)
    - +, unary - and - have been overloaded to represent S1 group operations, and always return S1 objects. 
    - Casts naturally to double
    - Already existing S1 objects can be assigned to a double (es. S1 a; a =
      50.2;), while construction from doubles must be explicit


  - local.cpp: Free particle on a 1D circle with simple local Metropolis
algorithm.   
  This algorithm is not ideal, as it presents an exponential critical slowing down for η -> 0 (issues with ergodicity)
  - tailor.cpp: Free particle on a 1D circle with local Metropolis + nonlocal
  "Tailor" move  
  Reference: C. Bonati, M. D'Elia, _Topological critical slowing down:
    Variations on a toy model_, Phys. Rev. E 98, 013308 (2018)
  - run.sh: Simple launcher script with GNU Parallel
