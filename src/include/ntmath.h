/**
 * @file ntmath.h
 * @author Maxim Deryabin (maxim.deryabin@gmail.com)
 * @brief Low level functions for math and Number Theory.
 * @version 0.1
 * @date 2023-05-14
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef __NTMATH_H__
#define __NTMATH_H__

#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <exception>

namespace polycalc {

typedef unsigned __int128 uint128_t;

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
inline void ModMultBarrettEq(uint64_t &a, uint64_t b, uint64_t m, uint64_t prec, size_t logm) {
    if (logm == 0) logm = MSB(m);
    if (prec == 0) prec = BarrettPrecompute(m, logm);

    uint128_t mul = static_cast<uint128_t>(a) * b;

    uint128_t tmp1 = mul;
    uint128_t tmp2 = tmp1 >> (logm-2);

    tmp1 = tmp2 * prec;
    tmp2 = tmp1 >> (logm + 5);
    tmp1 = tmp2 * m;

    a = static_cast<uint64_t>(mul - tmp1);

    while (a >= m) a -= m;
}

uint64_t ShoupPrecompute(uint64_t c, uint64_t m);
uint64_t ModMulShoup(uint64_t a, uint64_t c, uint64_t m, uint64_t prec);
inline void ModMulShoupEq(uint64_t &a, uint64_t c, uint64_t m, uint64_t prec) {
    // uint128_t aa  = static_cast<uint128_t>(a);
    uint128_t mul = static_cast<uint128_t>(a) * c;
    uint128_t tmp = ((static_cast<uint128_t>(a) * prec) >> 64) * m;

    a = static_cast<uint64_t>(mul - tmp);  // mod 2^64

    if (a >= m) a -= m;
}

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