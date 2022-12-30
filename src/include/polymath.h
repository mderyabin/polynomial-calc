#ifndef __POLYMATH_H__
#define __POLYMATH_H__

#include "math.h"

void ModAdd(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m);
void ModAdd(uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m);

void NaiveNegacyclicConvolution(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m);
void NaiveNegacyclicConvolution(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m, uint64_t mu, uint64_t logm);

#endif /* __POLYMATH_H__ */