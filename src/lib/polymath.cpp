#include <iostream>
#include <random>
#include <algorithm>

#include "polymath.h"

namespace polycalc {

void GenerateUniformPoly(uint64_t *ax, size_t N, uint64_t m) {
    std::random_device rd;
    std::default_random_engine en(rd());
    std::uniform_int_distribution<int> rnd(0, m-1);

    for (size_t i = 0; i < N; i++) { 
        ax[i] = rnd(en);
    }
}

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

void ModHadamardMul(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m, uint64_t prec, size_t logm) {
    if (logm == 0) logm = MSB(m);
    if (prec == 0) prec = BarrettPrecompute(m, logm);
    for (size_t i = 0; i < N; i++) {
        cx[i] = ModMultBarrett(ax[i], bx[i], m, prec, logm);
    }
}

void ModHadamardMul(uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m, uint64_t prec, size_t logm) {
    if (logm == 0) logm = MSB(m);
    if (prec == 0) prec = BarrettPrecompute(m, logm);
    for (size_t i = 0; i < N; i++) {
        ModMultBarrettEq(ax[i], bx[i], m, prec, logm);
    }
}

uint64_t ComputeValue(uint64_t x, const uint64_t *ax, size_t N, uint64_t m) {
    uint64_t res = ax[0];
    uint64_t xpow = x;
    uint64_t pr = ShoupPrecompute(x, m);
    size_t logm = MSB(m);
    uint64_t mu = BarrettPrecompute(m, logm);
    for (size_t i = 1; i < N; i++) {
        res = ModAdd(res, ModMultBarrett(ax[i], xpow, m, mu, logm), m);
        xpow = ModMulShoup(xpow, x, m, pr);
    }
    return res;
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

size_t ReverseBits(size_t number, size_t bit_length) {
    // Reverses the bits of `number` up to `bit_length`.
    size_t reversed = 0;
    for (size_t i = 0; i < bit_length; i++) {
        if ((number >> i) & 1)
            reversed |= 1 << (bit_length - 1 - i);
    }
    return reversed;
}

// twiddle factors generated here absorb powers of 2N-th roots to avoid NWC 
void ComputeTwiddleFactors(uint64_t *tf, size_t N, uint64_t m, bool isinverse) {
    uint64_t g = FindGenerator(m, N << 1);

    if (isinverse) 
        g = ModInvPrime(g, m);

    uint64_t *tf_direct = new uint64_t[N];
    size_t logN = MSB(N) - 1;

    uint64_t gg = ModMult(g, g, m);
    gg = ModMult(gg, g, m);

    uint64_t prec = ShoupPrecompute(gg, m);

    tf_direct[0] = 1;
    for (size_t i = 1; i < N; i++) {
        tf_direct[i] = ModMulShoup(tf_direct[i-1], gg, m, prec);
    }

    for (size_t i = 0; i < N; i++) {
        tf[i] = tf_direct[ReverseBits(i, logN)];
    }
    
    delete [] tf_direct;
}

void ComputeTwiddleFactorsNaive(uint64_t *tf, size_t N, uint64_t m, bool isinverse) {
    uint64_t g = FindGenerator(m, N << 1);

    if (isinverse) 
        g = ModInvPrime(g, m);

    uint64_t gg = ModMult(g, g, m);
    uint64_t prec = ShoupPrecompute(gg, m);

    tf[0] = 1;
    for (size_t i = 1; i < N; i++) {
        tf[i] = ModMulShoup(tf[i-1], gg, m, prec);
    }
}

void ComputeNWCSequence(uint64_t *pows, size_t N, uint64_t m, bool isinverse) {
    uint64_t g = FindGenerator(m, 2*N);
    if (isinverse) g = ModInvPrime(g, m);

    uint64_t prec = ShoupPrecompute(g, m);

    pows[0] = 1;
    for (size_t i = 1; i < N; i++) {
        pows[i] = ModMulShoup(pows[i-1], g, m, prec);
    }
}

void NWC(uint64_t *res, const uint64_t *ax, const uint64_t *pows, size_t N, uint64_t m) {
    size_t logm = MSB(m);
    uint64_t prec = BarrettPrecompute(m, logm);

    for (size_t i = 1; i < N; i++) {
        res[i] = ModMultBarrett(ax[i], pows[i], m, prec, logm);
    }
}

void iNWC(uint64_t *res, const uint64_t *ax, const uint64_t *ipows, size_t N, uint64_t m) {
    size_t logm = MSB(m);
    uint64_t prec = BarrettPrecompute(m, logm);

    for (size_t i = 1; i < N; i++) {
        res[i] = ModMultBarrett(ax[i], ipows[i], m, prec, logm);
    }
}

void NaiveNTT(uint64_t *res, const uint64_t *ax, const uint64_t *tf, const uint64_t *pows, size_t N, uint64_t m) {
    uint64_t *tx = new uint64_t[N];
    NWC(tx, ax, pows, N, m);
    for (size_t i = 0; i < N; i++) {
        res[i] = ComputeValue(tf[i], tx, N, m);
    }
    delete [] tx;
}

void NaiveInvNTT(uint64_t *res, const uint64_t *ax, const uint64_t *itf, const uint64_t *ipows, size_t N, uint64_t m) {
    uint64_t *tx = new uint64_t[N];
    uint64_t invN = ModInvPrime(N, m);
    uint64_t prec = ShoupPrecompute(invN, m);
    for (size_t i = 0; i < N; i++) {
        tx[i] = ComputeValue(itf[i], ax, N, m);
        ModMulShoupEq(tx[i], invN, m, prec);
    }
    iNWC(res, tx, ipows, N, m);
    delete [] tx;
}

// using algorithm from https://eprint.iacr.org/2016/504.pdf
void CooleyTukeyForwardNTT(uint64_t *ax, const uint64_t *tf, size_t N, uint64_t m, uint64_t prec, size_t logm) {
    size_t t = N;
    for (size_t n = 1; n < N; n <<= 1) {
        t >>= 1;
        for (size_t i = 0; i < n; i++) {
            size_t j1 = 2 * i * t;
            size_t j2 = j1 + t - 1;
            uint64_t s = tf[n + i]; 
            for (size_t j = j1; j <= j2; j++) {
                uint64_t u = ax[j];
                uint64_t v = ModMultBarrett(ax[j + t], s, m, prec, logm);
                ax[j] = ModAdd(u, v, m);
                ax[j + t] = ModSub(u, v, m);
            }   
        }
    }
}

// using algorithm from https://eprint.iacr.org/2016/504.pdf
void GentlemanSandeInverseNTT(uint64_t *ax, const uint64_t *itf, size_t N, uint64_t m, uint64_t invN, uint64_t prec_b, uint64_t prec_s, size_t logm) {
    size_t t = 1;
    for (size_t h = N>>1; h > 0; h >>= 1) {
        size_t j1 = 0;
        for (size_t i = 0; i < h; i++) {
            size_t j2 = j1 + t - 1;
            uint64_t s = itf[h + i]; 
            for (size_t j = j1; j <= j2; j++) {
                uint64_t u = ax[j];
                uint64_t v = ax[j + t];
                ax[j] = ModAdd(u, v, m);
                ax[j + t] = ModSub(u, v, m);
                ModMultBarrettEq(ax[j + t], s, m, prec_b, logm);
            }
            j1 += (t<<1);
        }
        t <<= 1;
    }

    for (size_t j = 0; j < N; j++) {
        ModMulShoupEq(ax[j], invN, m, prec_s);
    }
}

}