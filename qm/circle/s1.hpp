#ifndef S1_H
#define S1_H

#include <iostream>
class S1 {
    private:
        double x;
        void wrap();
    public:
        S1() :x{0} {}
        S1(const S1&);
        explicit S1(double a) { x = a; wrap(); }
        operator double() const { return x; }
        S1& operator =(double);
        S1& operator =(S1);
        S1& operator +=(double);
        S1& operator -=(double);
        std::ostream& operator<< (std::ostream&);
};

S1 operator +(double, S1);
S1 operator +(S1, double);
S1 operator +(S1, S1);
S1 operator -(S1, double);
S1 operator -(S1);


#endif
