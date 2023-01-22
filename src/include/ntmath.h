#ifndef __NTMATH_H__
#define __NTMATH_H__

#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <exception>

namespace polycalc {

size_t MSB(uint64_t t);

// return a + b mod m
uint64_t ModAdd(const uint64_t a, const uint64_t b, const uint64_t m);
void ModAddEq(uint64_t &a, const uint64_t b, const uint64_t m);

//
uint64_t ModSub(const uint64_t a, const uint64_t b, const uint64_t m);
void ModSubEq(uint64_t &a, const uint64_t b, const uint64_t m);

uint64_t ModMult(uint64_t a, uint64_t b, uint64_t m);
void ModMultEq(uint64_t &a, uint64_t b, uint64_t m);

// return a * b mod
uint64_t BarrettPrecompute(uint64_t m, size_t logm);
uint64_t ModMultBarrett(uint64_t a, uint64_t b, uint64_t m, uint64_t prec, size_t logm);
void ModMultBarrettEq(uint64_t &a, uint64_t b, uint64_t m, uint64_t prec, size_t logm);

uint64_t ShoupPrecompute(uint64_t c, uint64_t m);
uint64_t ModMulShoup(uint64_t a, uint64_t c, uint64_t m, uint64_t prec);
void ModMulShoupEq(uint64_t &a, uint64_t c, uint64_t m, uint64_t prec);

uint64_t ModExp(uint64_t a, uint64_t e, uint64_t m, uint64_t prec = 0, size_t logm = 0);

bool IsPrime(uint64_t m, size_t iters = 1000); // Miller-Rabin primarity test

uint64_t FindFirstPrimeDown(size_t logm, size_t M); // find prime in form m = (2^logm + 1) - k*M
uint64_t FindFirstPrimeUp(size_t logm, size_t M); // find prime in form m = (2^logm + 1) - k*M

uint64_t FindPrevPrime(uint64_t m, size_t M); // find previous prime in form m1 = m - k*M
uint64_t FindNextPrime(uint64_t m, size_t M); // find next prime in form m = m + k*M

uint64_t GCD(uint64_t a, uint64_t b);

uint64_t RhoPollard(uint64_t a);

std::vector<uint64_t> Factorize(uint64_t a);

uint64_t FindPrimitive(uint64_t n);
uint64_t FindGenerator(uint64_t m, size_t M);

uint64_t ModInvPrime(uint64_t a, uint64_t m);

}

#endif /* __NTMATH_H__ */