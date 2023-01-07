#ifndef __ENCODE_H__
#define __ENCODE_H__

#include "polynomial.h"

#include <cstdlib>
#include <vector>
#include <complex>

typedef std::complex<long double> cdouble;
typedef std::vector<cdouble> cvector;


Polynomial encode(cvector values, size_t N, uint64_t m, uint64_t scale); 

cvector decode(Polynomial poly, uint64_t scale);

#endif /* __ENCODE_H__ */