/**
 * @file polymath.h
 * @author Maxim Deryabin (maxim.deryabin@gmail.com)
 * @brief Some low level polynomial functions, including NTT/INTT
 * @version 0.1
 * @date 2023-05-14
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef __POLYMATH_H__
#define __POLYMATH_H__

#include "ntmath.h"

namespace polycalc {

void GenerateUniformPoly(uint64_t *ax, size_t N, uint64_t m);
void GenerateDiscreteGaussPoly(uint64_t *ax, size_t N, uint64_t m, double sigma = 3.19);
void GenerateBinaryPoly(uint64_t *ax, size_t N);

void ModNegate(uint64_t *ax, size_t N, uint64_t m);
void ModNegate(uint64_t *cx, uint64_t *ax, size_t N, uint64_t m);

void ModAdd(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m);
void ModAdd(uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m);

void ModSub(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m);
void ModSub(uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m);

void ModMul(uint64_t *cx, const uint64_t *ax, const uint64_t b, size_t N, uint64_t m);
void ModMul(uint64_t *ax, const uint64_t b, size_t N, uint64_t m);

void ModMul(uint64_t *cx, const uint64_t *ax, const uint64_t b, size_t N, uint64_t m, uint64_t mu, size_t logm);
void ModMul(uint64_t *ax, const uint64_t b, size_t N, uint64_t m, uint64_t mu, size_t logm);

void ModHadamardMul(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m, uint64_t prec = 0, size_t logm = 0);
void ModHadamardMul(uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m, uint64_t prec = 0, size_t logm = 0);

uint64_t ComputeValue(uint64_t x, const uint64_t *ax, size_t N, uint64_t m);

void NaiveNegacyclicConvolution(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m);
void NaiveNegacyclicConvolution(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m, uint64_t mu, uint64_t logm);

size_t ReverseBits(size_t number, size_t bit_length);
void ComputeTwiddleFactors(uint64_t *tf, size_t N, uint64_t m, bool isinverse = false);
void ComputeTwiddleFactorsNaive(uint64_t *tf, size_t N, uint64_t m, bool isinverse = false);
void ComputeNWCSequence(uint64_t *pows, size_t N, uint64_t m, bool isinverse = false);

void NaiveNTT(uint64_t *res, const uint64_t *ax, const uint64_t *tf, const uint64_t *pows, size_t N, uint64_t m);
void NaiveInvNTT(uint64_t *res, const uint64_t *ax, const uint64_t *itf, const uint64_t *ipows, size_t N, uint64_t m);

void CooleyTukeyForwardNTT(uint64_t *ax, const uint64_t *tf, size_t N, uint64_t m, uint64_t prec, size_t logm);
void GentlemanSandeInverseNTT(uint64_t *ax, const uint64_t *itf, size_t N, uint64_t m, uint64_t invN, uint64_t prec_b, uint64_t prec_s, size_t logm);

void ShoupPrecompute(uint64_t *prec_c, const uint64_t *c, size_t n, uint64_t m);
void CooleyTukeyForwardNTT(uint64_t *ax, const uint64_t *tf, size_t N, uint64_t m, const uint64_t *prec_tf, size_t logN);
void GentlemanSandeInverseNTT(uint64_t *ax, const uint64_t *itf, size_t N, uint64_t m, uint64_t invN, const uint64_t *prec_itf, uint64_t prec_invN);

void StockhamNTT(uint64_t *ax, const uint64_t *tf, size_t N, uint64_t q, uint64_t prec, size_t logm);
void ComputeTwiddleStockham(uint64_t *tf, size_t N, uint64_t m, bool isinverse);
void NWC(uint64_t *res, const uint64_t *ax, const uint64_t *pows, size_t N, uint64_t m);
} 

#endif /* __POLYMATH_H__ */