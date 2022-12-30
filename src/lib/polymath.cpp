#include "polymath.h"

void ModAdd(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m) {
    for (size_t i = 0; i < N; i++) {
        cx[i] = ModAdd(ax[i], bx[i], m);
    }
}

void ModAdd(uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m) {
    for (size_t i = 0; i < N; i++) {
        ModAddEq(ax[i], bx[i], m);
    }
}

void NaiveNegacyclicConvolution(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m) {
    uint64_t logm = MSB(m);
    uint64_t mu = BarrettPrecompute(m, logm);
    NaiveNegacyclicConvolution(cx, ax, bx, N, m, mu, logm);
}

void NaiveNegacyclicConvolution(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m, uint64_t mu, uint64_t logm) {
    for (size_t i = 0; i < N; i++) {
        cx[i] = 0;
        for (size_t j = 0; j <= i; j++) {
            // cx[i] += ax[j] * bx[i - j];
            ModAddEq(cx[i], ModMultBarrett(ax[j], bx[i - j], m, mu, logm), m);
        }
        for (size_t j = i+1; j <= N-1; j++) {
            // cx[i] -= ax[j] * bx[N + i - j];
            ModSubEq(cx[i], ModMultBarrett(ax[j], bx[N + i - j], m, mu, logm), m);
        }
    }   
}