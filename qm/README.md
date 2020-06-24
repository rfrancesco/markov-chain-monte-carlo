# QM simulations based on Metropolis algorithm


- harmonic.c: Quantum harmonic oscillator.
    
   
- circle/: Quantum free particle on a 1D circle  
  Contents:
  - pathsim.hpp/cpp: C++ class for simulating a free particle on a 1D circle with MCMC
    - The engine of these simulations
    - Measurements implemented as public functions
    - MCMC algorithms implemented as public functions
    - With some minor tweaks to the code I will eventually make, it is easily composable
      to implement parallel algorithms such as Parallel Tempering
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
  - slab.cpp: Free particle on a 1D circle. Technique for extracting the topological susceptibility with a local algorithm,
  when ergodicity among different topological sectors is lost (topological charge is frozen at Q = 0).
  - run.sh: Simple launcher script with GNU Parallel
