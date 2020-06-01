#ifndef S1_H
#define S1_H

#include <iostream>
class S1 {
    private:
        double x;
        void wrap();
    public:
	// Constructor / Destructor
        S1() :x{0} {}
        S1(const S1&);
        explicit S1(double a) { x = a; wrap(); }
	// Cast to double
        operator double() const { return x; }
	// Assignment and group operations (+, -)
        S1& operator =(double);
        S1& operator =(S1);
        S1& operator +=(double);
        S1& operator -=(double);
	// Print to ostream
        std::ostream& operator<< (std::ostream&);
};

// Group operations (+, -) with S1 and with doubles
S1 operator +(double, S1);
S1 operator +(S1, double);
S1 operator +(S1, S1);
S1 operator -(S1, double);
S1 operator -(S1);


#endif
