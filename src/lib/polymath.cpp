#include "polymath.h"

void ModAdd(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m) {
    for (size_t i = 0; i < N; i++) {
        cx[i] = ModAdd(ax[i], bx[i], m);
    }
}

void ModAdd(const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m) {
    for (size_t i = 0; i < N; i++) {
        ModAddEq(ax[i], bx[i], m);
    }
}

void NaiveNegacyclicConvolution(int64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m) {
    // TODO: реализовать умножение
}