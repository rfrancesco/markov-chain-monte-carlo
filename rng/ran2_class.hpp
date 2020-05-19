#ifndef RAN2_H
#define RAN2_H

#include <cfloat>
#include <cmath>
#include <iostream>

// Ran2 class header
// This is a dummy header without the accompanying .cpp file.

#define NTAB 32

class Ran2 {
    private:
    // Internal state
        std::int64_t idum;
        std::int64_t idum2;
        std::int64_t iy;
        std::int64_t iv[NTAB];
    public:
        // Constructor - takes integer seed > 0 (sign is flipped inside the constructor)
        explicit Ran2(std::int64_t);
        // Generate random double (uniform deviates)
        double doub();
        // Generate gaussian double with Box-Muller algorithm
        double gaussian(double, double);
};

#undef NTAB


#endif
