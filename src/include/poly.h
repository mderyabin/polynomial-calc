#ifndef __POLY_H__
#define __POLY_H__

#include <cstdint>
#include <cstddef>
#include <iostream>

#include "math.h"

class Poly {
private:
    uint64_t *ax;  // coefficients
    uint64_t m;    // modulus
    size_t N;      // degree

    uint64_t mu;   // Barrett mul precompute constant 
    uint64_t logm; // size of modulus in bits
public:
    // constructors
    Poly() : ax(NULL), N(0), m(0) { }
    Poly(size_t _N, uint64_t _m, bool initWithZeros = true);
    Poly(uint64_t *_ax, size_t _N, uint64_t _m);

    Poly(const Poly &o); // copy
    Poly& operator=(const Poly &o);

    Poly(Poly &&o); // move
    Poly& operator=(Poly &&o);
    
    // destructor
    ~Poly() { if (ax) delete [] ax; }

    // get parameters
    size_t GetN() const { return N; }
    uint64_t GetModulus() const { return m; }

    // indexing 
    uint64_t &operator[](size_t i) { return ax[i]; }
    uint64_t at(size_t i) const {return ax[i]; }

    // calculate value of the polynomial
    uint64_t operator()(uint64_t x) const;

    // print to the stream
    friend std::ostream& operator<<(std::ostream& os, const Poly& poly);

    // arithmetic 
    friend const Poly operator+(const Poly& left, const Poly& right);
    friend const Poly& operator+=(Poly& left, const Poly& right);
    friend const Poly& operator+(const Poly& poly) {return poly; }
};

#endif /* __POLY_H__ */