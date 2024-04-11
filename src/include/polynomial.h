/**
 * @file polynomial.h
 * @author Maxim Deryabin (maxim.deryabin@gmail.com)
 * @brief Class for organization computations with (non-rns) native polynomial.
 * @version 0.1
 * @date 2023-05-14
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include <iostream>
#include <memory>
#include <string>

#include "polymath.h"
#include "utils.h"
#include "ntt.h"

#include "cereal.hpp"


namespace polycalc {

/**
 * @brief Determines the state of polynomial
 * 
 * EVAL if polynomial is in NTT format
 * COEF if polynomial is represented in coefficients
 */
enum Format { EVAL, COEF };

/**
 * @brief Polynomial class which supports arithmetic polynomial operations.
 * 
 * Internally polynomial is stored as native C array. 
 * This class includes NTT/INTT. 
 */
class Polynomial {
private:
    uint64_t *ax = NULL;     // coefficients
    size_t N;                // degree
    uint64_t m;              // modulus
    uint64_t logm;           // size of modulus in bits
    uint64_t mu;             // Barrett mul precompute constant 

    Format format;

    std::shared_ptr<NTT> ntt;

    friend class cereal::access;

    // order: format (0 for EVAL and 1 for COEF), N, m, logm, mu, ax[i]
    template<class Archive>
    void save(Archive & archive) const;

    template<class Archive>
    void load(Archive & archive);
public:
    /****** constructors and destructors and assigment operators ******/
    /******************************************************************/
    Polynomial() : ax(NULL), N(0), m(0), format(COEF) { }
    Polynomial(size_t _N, uint64_t _m, bool initWithZeros = true, Format _format = COEF);
    Polynomial(uint64_t *_ax, size_t _N, uint64_t _m, Format _format = COEF);

    Polynomial(const Polynomial &o); // copy
    Polynomial& operator=(const Polynomial &o);

    Polynomial(Polynomial &&o); // move
    Polynomial& operator=(Polynomial &&o);
    
    ~Polynomial() { if (ax) delete [] ax; }
    /******************************************************************/

    /****** general polynomial operations ******/
    // Fill polynomial with random values
    void GenerateUniform(Format _format);
    void GenerateUniform();

    // fill with Discrete Gauss
    // change format to coeff 
    void GenerateDiscreteGauss();

    void GenerateBinary();

    // Fill polynomial with zeroes
    void SetZero(Format _format);
    void SetZero();

    // Calculate value of the Polynomial on value x
    uint64_t operator()(uint64_t x) const;
    /*******************************************/

    /****** access to parameters ******/
    size_t GetN() const { return N; }
    uint64_t GetModulus() const { return m; }
    Format GetFormat() const { return format; }
    /**********************************/

    /****** formats ******/
    void SetFormatEval();
    void SetFormatCoef();
    /*********************/

    /****** arithmetics ******/   
    void NegateInPlace();

    friend const Polynomial operator+(const Polynomial& left, const Polynomial& right);
    friend const Polynomial& operator+=(Polynomial& left, const Polynomial& right);
    friend const Polynomial& operator+(const Polynomial& Polynomial) {return Polynomial; }

    friend const Polynomial operator*(const Polynomial& left, const Polynomial& right);
    friend const Polynomial& operator*=(Polynomial& left, const Polynomial& right);
    /*************************/

    /****** input and output ******/
    // Save polynomial into file using serialization with cereal
    void Serialize(std::string filename, SER_Archive_Type TYPE = BIN);

    // Load polynomial from file using serialization with cereal
    void Deserialize(std::string filename, SER_Archive_Type TYPE = BIN);

    // Print to the stream in form of natural polynomial
    friend std::ostream& operator<<(std::ostream& os, const Polynomial& Polynomial);

    // Indexing 
    uint64_t &operator[](size_t i) { return ax[i]; }

    // Get coefficient/value by index
    uint64_t at(size_t i) const {return ax[i]; }
    /******************************/
};

}

#endif /* __POLYNOMIAL_H__ */