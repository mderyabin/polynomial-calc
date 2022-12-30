#ifndef __MATH_H__
#define __MATH_H__

#include <cstdint>
#include <cstddef>

uint64_t MSB(uint64_t t);

// return a + b mod m
uint64_t ModAdd(const uint64_t a, const uint64_t b, const uint64_t m);
void ModAddEq(uint64_t &a, const uint64_t b, const uint64_t m);

// return a * b mod
uint64_t BarrettPrecompute(uint64_t m, uint64_t logm);
uint64_t ModMultBarrett(uint64_t a, uint64_t b, uint64_t m, uint64_t prec, uint64_t logm);

uint64_t ShoupPrecompute(uint64_t c, uint64_t m);
uint64_t ModMulShoup(uint64_t a, uint64_t c, uint64_t m, uint64_t prec);



#endif /* __MATH_H__ */