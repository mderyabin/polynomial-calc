#include "benchmark/benchmark.h"

#include "polymath.h"

using namespace polycalc;

const size_t logN = 12;
const size_t logq = 58;

void BC_NaiveNTT(benchmark::State& state) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);


    uint64_t *tf = new uint64_t[N];
    uint64_t *pows = new uint64_t[N];
    uint64_t *ax = new uint64_t[N];
    uint64_t *nttax = new uint64_t[N];;

    GenerateUniformPoly(ax, N, q); 

    ComputeTwiddleFactorsNaive(tf, N, q, false);
    ComputeNWCSequence(pows, N, q, false);

    while (state.KeepRunning()) {
        NaiveNTT(nttax, ax, tf, pows, N, q);
    }

    delete [] nttax;
    delete [] ax;
    delete [] pows;
    delete [] tf;
}

BENCHMARK(BC_NaiveNTT)->Unit(benchmark::kNanosecond);

void BC_CooleyTukeyForwardNTT(benchmark::State& state) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);
    uint64_t prec = BarrettPrecompute(q, logq);


    uint64_t *tf = new uint64_t[N];
    uint64_t *ax = new uint64_t[N];

    GenerateUniformPoly(ax, N, q); 

    ComputeTwiddleFactorsNaive(tf, N, q, false);

    while (state.KeepRunning()) {
        CooleyTukeyForwardNTT(ax, tf, N, q, prec, logq);
    }

    delete [] ax;
    delete [] tf;
}

BENCHMARK(BC_CooleyTukeyForwardNTT)->Unit(benchmark::kNanosecond);

void BC_CooleyTukeyForwardNTT_Shoup(benchmark::State& state) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);

    uint64_t *tf = new uint64_t[N];
    uint64_t *prec_tf = new uint64_t[N];
    uint64_t *ax = new uint64_t[N];

    GenerateUniformPoly(ax, N, q); 
    ComputeTwiddleFactorsNaive(tf, N, q, false);
    ShoupPrecompute(prec_tf, tf, N, q);

    while (state.KeepRunning()) {
        CooleyTukeyForwardNTT(ax, tf, N, q, prec_tf, logN);
    }

    delete [] ax;
    delete [] prec_tf;
    delete [] tf;
}

BENCHMARK(BC_CooleyTukeyForwardNTT_Shoup)->Unit(benchmark::kNanosecond);

void BC_NaiveInvNTT(benchmark::State& state) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);


    uint64_t *itf = new uint64_t[N];
    uint64_t *ipows = new uint64_t[N];
    uint64_t *ax = new uint64_t[N];
    uint64_t *nttax = new uint64_t[N];;

    GenerateUniformPoly(ax, N, q); 

    ComputeTwiddleFactorsNaive(itf, N, q, true);
    ComputeNWCSequence(ipows, N, q, true);

    while (state.KeepRunning()) {
        NaiveInvNTT(nttax, ax, itf, ipows, N, q);
    }

    delete [] nttax;
    delete [] ax;
    delete [] ipows;
    delete [] itf;
}

BENCHMARK(BC_NaiveInvNTT)->Unit(benchmark::kNanosecond);

void BC_GentlemanSandeInverseNTT(benchmark::State& state) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);
    uint64_t prec = BarrettPrecompute(q, logq);

    uint64_t invN = ModInvPrime(N, q);
    uint64_t prec_invN = ShoupPrecompute(invN, q);


    uint64_t *itf = new uint64_t[N];
    uint64_t *ax = new uint64_t[N];

    GenerateUniformPoly(ax, N, q); 

    ComputeTwiddleFactorsNaive(itf, N, q, true);

    while (state.KeepRunning()) {
        GentlemanSandeInverseNTT(ax, itf, N, q, invN, prec, prec_invN, logq);
    }

    delete [] ax;
    delete [] itf;
}

BENCHMARK(BC_GentlemanSandeInverseNTT)->Unit(benchmark::kNanosecond);

void BC_GentlemanSandeInverseNTT_Shoup(benchmark::State& state) {
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);

    uint64_t invN = ModInvPrime(N, q);
    uint64_t prec_invN = ShoupPrecompute(invN, q);

    uint64_t *itf = new uint64_t[N];
    uint64_t *prec_itf = new uint64_t[N];
    uint64_t *ax = new uint64_t[N];

    GenerateUniformPoly(ax, N, q); 
    ComputeTwiddleFactorsNaive(itf, N, q, true);
    ShoupPrecompute(prec_itf, itf, N, q);

    while (state.KeepRunning()) {
        GentlemanSandeInverseNTT(ax, itf, N, q, invN, prec_itf, prec_invN);
    }

    delete [] ax;
    delete [] prec_itf;
    delete [] itf;
}

BENCHMARK(BC_GentlemanSandeInverseNTT_Shoup)->Unit(benchmark::kNanosecond);


BENCHMARK_MAIN();