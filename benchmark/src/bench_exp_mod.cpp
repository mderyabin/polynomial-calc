#include <iostream>


#include <NTL/ZZ.h>

#include "benchmark/benchmark.h"

#include "vectbigint.h"

using namespace std;
using namespace polycalc;
using namespace NTL;

const uint64_t q = (1ULL << 30) - (1ULL << 18) + 1;
const size_t s = 8;

void TransformZZtoCPoly(uint64_t *ax, const NTL::ZZ a, size_t n, size_t s) {    
    NTL::ZZ a_copy = a;
    for (size_t i = 0; i < n; i++) {
        ax[i] = a_copy % (1<<s);
        a_copy /= (1<<s);
    }
}

ZZ FastModExp(ZZ a, ZZ e, ZZ m) { 
    ZZ res, base, ee;
    bool flag = false;
    base = a;
    ee = e;

    while (ee > 0) {
        if (ee % 2 == 1) {
            if (flag) { 
                res = MulMod(res, base, m);
            } else {
                res = base;
                flag = true;
            }
        }
        ee >>= 1;
        if (ee != 0) {
            base = SqrMod(base, m);
        }
    }

    return res;
}

void BC_ZZ(benchmark::State& state, size_t N, size_t exp_size) {
    size_t bsize = s*N / 2;
    ZZ a, e, m, c;
    m = GenPrime_ZZ(bsize-2);
    a = RandomBits_ZZ(bsize) % m;
    e = RandomBits_ZZ(exp_size);

    while (state.KeepRunning()) {
        PowerMod(c, a, e, m);
    }
}

void BC_Schoolbook(benchmark::State& state, size_t N, size_t exp_size) {
    size_t bsize = s*N / 2;
    ZZ a, e, m;
    m = GenPrime_ZZ(bsize-2);
    a = RandomBits_ZZ(bsize) % m;
    e = RandomBits_ZZ(exp_size);


    uint64_t mx[N];

    uint64_t ax[N];
    uint64_t ex[N];
    uint64_t cx[N];


    TransformZZtoCPoly(ax, a, N, s);
    TransformZZtoCPoly(ex, e, N, s);
    TransformZZtoCPoly(mx, m, N, s);


    while (state.KeepRunning()) {
        FastModExp(cx, ax, ex, mx, N, s);
    }
}

void BC_REDC(benchmark::State& state, size_t N, size_t exp_size) {
    size_t bsize = s*N / 2;
    ZZ a, e, m;
    m = GenPrime_ZZ(bsize-2);
    a = RandomBits_ZZ(bsize) % m;
    e = RandomBits_ZZ(exp_size);

    size_t r_bits = bsize;
    uint64_t mu;

    uint64_t mx[N];
    uint64_t rsx[N];

    uint64_t ax[N];
    uint64_t ex[N];
    uint64_t cx[N];

    NTTInstance ntt = NTTManager::GetNTTPtr(N, q);
    REDC_precompute(mx, rsx, mu, m, r_bits, s, N);

    TransformZZtoCPoly(ax, a, N, s);
    TransformZZtoCPoly(ex, e, N, s);

    while (state.KeepRunning()) {
        FastModExpREDC(cx, ax, ex, mx, rsx, mu, r_bits, N, s, ntt);
    }
}

BENCHMARK_CAPTURE(BC_ZZ, BC_ZZ_128_128, 128, 128)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_ZZ, BC_ZZ_256_256, 256, 256)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_ZZ, BC_ZZ_512_512, 512, 512)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_ZZ, BC_ZZ_1024_1024, 1024, 1024)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_ZZ, BC_ZZ_2048_2048, 2048, 2048)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(BC_Schoolbook, BC_Schoolbook_128_128, 128, 128)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_Schoolbook, BC_Schoolbook_256_256, 256, 256)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_Schoolbook, BC_Schoolbook_512_512, 512, 512)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_Schoolbook, BC_Schoolbook_1024_1024, 1024, 1024)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_Schoolbook, BC_Schoolbook_2048_2048, 2048, 2048)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(BC_REDC, BC_REDC_128_128, 128, 128)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_REDC, BC_REDC_256_256, 256, 256)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_REDC, BC_REDC_512_512, 512, 512)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_REDC, BC_REDC_1024_1024, 1024, 1024)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BC_REDC, BC_REDC_2048_2048, 2048, 2048)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
