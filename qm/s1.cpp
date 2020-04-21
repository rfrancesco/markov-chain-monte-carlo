#include "s1.hpp"

S1::S1(const S1& value) {
    x = static_cast<double>(value);
}


void S1::wrap() {
    while (x < -0.5)
        x += 1;
    while (x >= 0.5)
        x -= 1;
}

S1& S1::operator =(double value) {
    x = value;
    wrap();
    return *this;
}

S1& S1::operator =(S1 value) {
    x = static_cast<double>(value);
    return *this;
}

S1& S1::operator +=(double value) {
    x += value;
    wrap();
    return *this;
}

S1& S1::operator -=(double value) {
    x -= value;
    wrap();
    return *this;
}

S1 operator +(double a, S1 b) {
    return b += a;
}

S1 operator +(S1 b, double a) {
    return b += a;
}

S1 operator +(S1 a, S1 b) {
    return b += a;
}

S1 operator -(S1 a, double b) {
    return a -= b;
}

S1 operator -(S1 a) {
    a = -static_cast<double>(a);
    return a;
}

std::ostream& S1::operator<< (std::ostream& stm) {
    return stm << x;
}
