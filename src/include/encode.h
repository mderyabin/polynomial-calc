/**
 * @file encode.h
 * @author Maxim Deryabin (maxim.deryabin@gmail.com)
 * @brief Functionality for encoding vector of values into Polynomial object in CKKS style. 
 * @version 0.1
 * @date 2023-05-14
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef __ENCODE_H__
#define __ENCODE_H__

#include "polynomial.h"

#include <cstdlib>
#include <vector>
#include <complex>

namespace polycalc {

typedef std::complex<long double> cdouble;
typedef std::vector<cdouble> cvector;

/**
 * @brief Encode vector of complex numbers into ring polynomial. 
 * 
 * This method uses special case of interpolation of complex-valued polynomial
 * which results in the polynomial with only real values. 
 * Using scaling and rounding, coefficients of that polynomial is reduced to integers modulo m.
 * 
 * @param values vector of complex values of the size <= N/2
 * @param N ring dimension
 * @param m modulus
 * @param scale smaller than m, used to scale and round real coefficients to integers
 * @return encoded ring Polynomial modulo m and degree N-1
 */
Polynomial encode(cvector values, size_t N, uint64_t m, uint64_t scale); 

cvector decode(Polynomial poly, uint64_t scale);

} 

#endif /* __ENCODE_H__ */