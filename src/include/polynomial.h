#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include <iostream>
#include <memory>

#include "polymath.h"
#include "ntt.h"

// determines the current state of polynomial
// EVAL if polynomial is in NTT format
// COEF if polynomial is represented in coefficients
enum Format { EVAL, COEF };

class Polynomial {
private:
    uint64_t *ax;  // coefficients
    uint64_t m;    // modulus
    size_t N;      // degree

    uint64_t mu;   // Barrett mul precompute constant 
    uint64_t logm; // size of modulus in bits

    Format format;

    std::shared_ptr<NTT> ntt;
public:
    // constructors
    Polynomial() : ax(NULL), N(0), m(0), format(COEF) { }
    Polynomial(size_t _N, uint64_t _m, bool initWithZeros = true, Format _format = COEF);
    Polynomial(uint64_t *_ax, size_t _N, uint64_t _m, Format _format = COEF);

    Polynomial(const Polynomial &o); // copy
    Polynomial& operator=(const Polynomial &o);

    Polynomial(Polynomial &&o); // move
    Polynomial& operator=(Polynomial &&o);
    
    // destructor
    ~Polynomial() { if (ax) delete [] ax; }

    // fill polynomial
    void GenerateUniform(Format _format = COEF);

    // get parameters
    size_t GetN() const { return N; }
    uint64_t GetModulus() const { return m; }
    Format GetFormat() const { return format; }

    // indexing 
    uint64_t &operator[](size_t i) { return ax[i]; }
    uint64_t at(size_t i) const {return ax[i]; }

    // calculate value of the Polynomialnomial
    uint64_t operator()(uint64_t x) const;

    // print to the stream
    friend std::ostream& operator<<(std::ostream& os, const Polynomial& Polynomial);

    // format change
    void SetFormatEval();
    void SetFormatCoef();

    // arithmetic 
    friend const Polynomial operator+(const Polynomial& left, const Polynomial& right);
    friend const Polynomial& operator+=(Polynomial& left, const Polynomial& right);
    friend const Polynomial& operator+(const Polynomial& Polynomial) {return Polynomial; }

    friend const Polynomial operator*(const Polynomial& left, const Polynomial& right);
    friend const Polynomial& operator*=(Polynomial& left, const Polynomial& right);
};

#endif /* __POLYNOMIAL_H__ */