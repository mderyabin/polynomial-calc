#include <iostream>
#include "encode.h"
#include "ntmath.h"
#include "polynomial.h"

#include "benchmark/benchmark.h"


#include <NTL/ZZ.h>

using namespace std;
using namespace polycalc;
using namespace NTL;

const size_t osize_bit = 1024;

const size_t N = 2*osize_bit;       // size of numbers (product of a and b size) 
const size_t w = 64;                // word size

const size_t M = 16;                 // coefficients size
const size_t k = N/M;               // number of coefficients

const size_t log_q = 58;          

void TransformZZtoCPoly(Polynomial &ax, const ZZ a) {
    ZZ a_copy = a;
    for (size_t i = 0; i < ax.GetN(); i++) {
        ax[i] = a_copy % (1<<M);
        a_copy /= (1<<M);
    }
}

void ReconstructZZfromCPoly(ZZ &a, const Polynomial ax) {
    a = 0;
    uint64_t x(1<<M);
    for (int i = ax.GetN() - 1; i >= 0; i--) {
        a = (a * x + ax.at(i));
    }
}

void BC_ZZBugNumMul(benchmark::State& state, size_t N) {
    ZZ a, b, c;
    a = RandomBits_ZZ(N/2);
    b = RandomBits_ZZ(N/2);

    while (state.KeepRunning()) {
        c = a * b;
    }
}

void BC_PolyBugNumMul(benchmark::State& state, size_t N, size_t M, size_t log_q) {
    size_t k = N/M;
    uint64_t q = FindFirstPrimeUp(log_q, 2*k);
    Polynomial ax(k, q);
    Polynomial bx(k, q);
    Polynomial cx(k, q);

    ZZ a, b, c;
    a = RandomBits_ZZ(N/2);
    b = RandomBits_ZZ(N/2);

    TransformZZtoCPoly(ax, a);
    TransformZZtoCPoly(bx, b);

    while (state.KeepRunning()) {
        ax.SetFormatEval();
        bx.SetFormatEval();
        cx = ax * bx;
        cx.SetFormatCoef();
    }
}

BENCHMARK_CAPTURE(BC_ZZBugNumMul, BC_ZZBugNumMul_256, 256);
BENCHMARK_CAPTURE(BC_ZZBugNumMul, BC_ZZBugNumMul_512, 512);
BENCHMARK_CAPTURE(BC_ZZBugNumMul, BC_ZZBugNumMul_1024, 1024);

BENCHMARK_CAPTURE(BC_PolyBugNumMul, BC_PolyBugNumMul_256_16_58, 256, 16, 58);
BENCHMARK_CAPTURE(BC_PolyBugNumMul, BC_PolyBugNumMul_512_16_58, 512, 16, 58);
BENCHMARK_CAPTURE(BC_PolyBugNumMul, BC_PolyBugNumMul_1024_16_58, 1024, 16, 58);

BENCHMARK_CAPTURE(BC_PolyBugNumMul, BC_PolyBugNumMul_256_8_58, 256, 8, 58);
BENCHMARK_CAPTURE(BC_PolyBugNumMul, BC_PolyBugNumMul_512_8_58, 512, 8, 58);
BENCHMARK_CAPTURE(BC_PolyBugNumMul, BC_PolyBugNumMul_1024_8_58, 1024, 8, 58);

BENCHMARK_MAIN();
