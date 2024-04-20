#include <random>

#include "benchmark/benchmark.h"

#include "polymath.h"

using namespace polycalc;

uint64_t generate_number(uint64_t m) {
    std::random_device rd;
    std::default_random_engine en(rd());
    std::uniform_int_distribution<int> rnd(0, m-1);

    return rnd(en);
}

void BC_ShoupMul(benchmark::State& state, size_t logN, size_t logq) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);

    uint64_t *ax = new uint64_t[N];
    GenerateUniformPoly(ax, N, q); 

    uint64_t b = generate_number(q);

    uint64_t *cx = new uint64_t[N];

    while (state.KeepRunning()) {
        ModMul(cx, ax, b, N, q);
    }

    delete [] ax;
    delete [] cx;
}

void BC_ShoupEqMul(benchmark::State& state, size_t logN, size_t logq) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);

    uint64_t *ax = new uint64_t[N];
    GenerateUniformPoly(ax, N, q); 

    uint64_t b = generate_number(q);

    while (state.KeepRunning()) {
        ModMul(ax, b, N, q);
    }

    delete [] ax;
}

void BC_BarrettMult(benchmark::State& state, size_t logN, size_t logq) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);

    uint64_t *ax = new uint64_t[N];
    GenerateUniformPoly(ax, N, q); 

    uint64_t b = generate_number(q);

    uint64_t *cx = new uint64_t[N];

    uint64_t mu = BarrettPrecompute(q, logq);

    while (state.KeepRunning()) {
        ModMul(cx, ax, b, N, q, mu, logq);
    }

    delete [] ax;
    delete [] cx;
}

void BC_BarrettEqMult(benchmark::State& state, size_t logN, size_t logq) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);

    uint64_t *ax = new uint64_t[N];
    GenerateUniformPoly(ax, N, q); 

    uint64_t b = generate_number(q);

    uint64_t mu = BarrettPrecompute(q, logq);

    while (state.KeepRunning()) {
        ModMul(ax, b, N, q, mu, logq);
    }

    delete [] ax;
}

BENCHMARK_CAPTURE(BC_ShoupMul, Shoup_12_58, 12, 58)->Unit(benchmark::kNanosecond);
BENCHMARK_CAPTURE(BC_ShoupEqMul, Shoup_12_58, 12, 58)->Unit(benchmark::kNanosecond);
BENCHMARK_CAPTURE(BC_BarrettMult, Barrett_12_58, 12, 58)->Unit(benchmark::kNanosecond);
BENCHMARK_CAPTURE(BC_BarrettEqMult, BarrettEq_12_58, 12, 58)->Unit(benchmark::kNanosecond);

BENCHMARK_MAIN();