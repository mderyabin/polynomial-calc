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
#define __ENCODE_H__

#include "polynomial.h"

#include <cstdlib>
#include <vector>
#include <complex>

namespace polycalc {

typedef std::complex<long double> cdouble;
typedef std::vector<cdouble> cvector;


Polynomial encode(cvector values, size_t N, uint64_t m, uint64_t scale); 

cvector decode(Polynomial poly, uint64_t scale);

} 

#endif /* __ENCODE_H__ */