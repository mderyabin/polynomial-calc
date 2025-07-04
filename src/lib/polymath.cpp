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

void GenerateBinaryPoly(uint64_t *ax, size_t N) {
    std::random_device rd;
    std::default_random_engine en(rd());
    std::uniform_int_distribution<int> rnd(0, 1);

    for (size_t i = 0; i < N; i++) { 
        ax[i] = rnd(en);
    }
}

void GenerateDiscreteGaussPoly(uint64_t *ax, size_t N, uint64_t m, double sigma) {
    std::random_device rd {};
    std::mt19937 gen {rd()};
    std::normal_distribution<> d {0, sigma};

    for (size_t i = 0; i < N; i++) { 
        int64_t t = std::round(d(gen));
        ax[i] = t >= 0 ? t : m + t;
    }
}

void ModNegate(uint64_t *ax, size_t N, uint64_t m) {
    for (size_t i = 0; i < N; i++) {
        ModNegateEq(ax[i], m);
    }
}

void ModNegate(uint64_t *cx, uint64_t *ax, size_t N, uint64_t m) {
    for (size_t i = 0; i < N; i++) {
        cx[i] = ModNegate(ax[i], m);
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

void ModSub(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m) {
    for (size_t i = 0; i < N; i++) {
        cx[i] = ModSub(ax[i], bx[i], m);
    }
}

void ModSub(uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m) {
    for (size_t i = 0; i < N; i++) {
        ModSubEq(ax[i], bx[i], m);
    }
}

void ModMul(uint64_t *cx, const uint64_t *ax, const uint64_t b, size_t N, uint64_t m) {
    uint64_t prec = ShoupPrecompute(b, m); 
    for (size_t i = 0; i < N; i++) {
        cx[i] = ModMulShoup(ax[i], b, m, prec);
    }
}

void ModMul(uint64_t *cx, const uint64_t *ax, const uint64_t b, size_t N, uint64_t m, uint64_t mu, size_t logm) {
    for (size_t i = 0; i < N; i++) {
        cx[i] = ModMultBarrett(ax[i], b, m, mu, logm);
    }
}

void ModMul(uint64_t *ax, const uint64_t b, size_t N, uint64_t m) {
    uint64_t prec = ShoupPrecompute(b, m); 
    for (size_t i = 0; i < N; i++) {
        ModMulShoupEq(ax[i], b, m, prec);
    }
}

void ModMul(uint64_t *ax, const uint64_t b, size_t N, uint64_t m, uint64_t mu, size_t logm) {
    for (size_t i = 0; i < N; i++) {
        ModMultBarrettEq(ax[i], b, m, mu, logm);
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

    // uint64_t gg = ModMult(g, g, m);
    // gg = ModMult(gg, g, m);

    uint64_t prec = ShoupPrecompute(g, m);

    tf_direct[0] = 1;
    for (size_t i = 1; i < N; i++) {
        tf_direct[i] = ModMulShoup(tf_direct[i-1], g, m, prec);
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

void ComputeTwiddleStockham(uint64_t *tf, size_t N, uint64_t m, bool isinverse) {
    uint64_t g = FindGenerator(m, N << 1);

    if (isinverse) 
        g = ModInvPrime(g, m);

    size_t logN = MSB(N) - 1;

    // uint64_t gg = ModMult(g, g, m);
    // gg = ModMult(gg, g, m);

    uint64_t prec = ShoupPrecompute(g, m);

    tf[0] = 1;
    for (size_t i = 1; i < N; i++) {
        tf[i] = ModMulShoup(tf[i-1], g, m, prec);
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
    // NWC(tx, ax, pows, N, m);
    for (size_t i = 0; i < N; i++) {
        // res[i] = ComputeValue(tf[i], tx, N, m);
        res[i] = ComputeValue(tf[i], ax, N, m);
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

void StockhamNTT(uint64_t *ax, const uint64_t *tf, size_t N, uint64_t q, uint64_t prec, size_t logq) { 
    uint64_t *buffer = new uint64_t[N];

    uint64_t *src, *dst;
    src = ax;
    dst = buffer;

    size_t hN = N>>1; // half of N

    size_t logN = MSB(N) - 1;

    size_t m, t, stride;
    size_t sidx1, sidx2;
    size_t didx1, didx2;

    uint64_t w, tw_v;
    uint64_t v, u;

    for (size_t s = 0; s < logN; s++) {
        m = 1<<s;               // Number of blocks
        t = N / (2 * m);        // Size of each butterfly group
        stride = N >> (s + 1);   // Distance between twiddles

        // std::cout << "s = " << s << std::endl; 
        // std::cout << "m = " << m << std::endl; 
        // std::cout << "t = " << t << std::endl; 
        // std::cout << "stride = " << stride << std::endl; 

        for (size_t j = 0; j < m; j++) {
            w = tf[j + stride]; // Twiddle for this butterfly

            // std::cout << " ----- " << std::endl;
            // std::cout << "w = " << w << std::endl << std::endl; 


            for (size_t i = 0; i < t; i++) {
                // compute indexes for source
                sidx1 = (2*j + 0)*t + i;
                sidx2 = (2*j + 1)*t + i;

                // compute indexes for destination
                didx1 = (j + 0) * t + i;
                didx2 = (j + m) * t + i;

                // std::cout << sidx1 << " " << sidx2 << std::endl; 
                // std::cout << didx1 << " " << didx2 << std::endl; 
                // std::cout << std::endl; 

                // read data
                u = src[sidx1];
                v = src[sidx2];

                // main mult
                ModMultBarrettEq(v, w, q, prec, logq);
                
                // butterfly add/sub
                dst[didx1] = ModAdd(u, v, q);
                dst[didx2] = ModSub(u, v, q);
            }
        }

        //swapping source and destination for the next stage
        std::swap(src, dst);
    }
    
    if (ax!=src)
        std::copy(src, src+N, ax);
    
    delete [] buffer;
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

void ShoupPrecompute(uint64_t *prec_c, const uint64_t *c, size_t n, uint64_t m) {
    for (size_t i = 0; i < n; i++) {
        prec_c[i] = ShoupPrecompute(c[i], m);
    }
    
}

void CooleyTukeyForwardNTT(uint64_t *ax, const uint64_t *tf, size_t N, uint64_t m, const uint64_t *prec_tf, size_t logN) {
    size_t t = N >> 1;
    size_t log1 = logN;
    uint64_t temp1, u, v, s, s_prec;
    size_t n, j1, j2, j, i;
    for (n = 1; n < N; n <<= 1) {
        for (i = 0; i < n; i++) {
            j1 = i << log1;
            j2 = j1 + t;
            s = tf[n + i]; 
            s_prec = prec_tf[n + i]; 
            for (j = j1; j < j2; j++) {
                u = ax[j];
                v = ax[j + t];

                ModMulShoupEq(v, s, m, s_prec);

                temp1 = u + v;
                ax[j] = temp1 < m ? temp1 : temp1 - m;

                ax[j + t] = (u > v) ? u - v : m + u - v;
            }   
        }
        t >>= 1;
        log1--;
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

void GentlemanSandeInverseNTT(uint64_t *ax, const uint64_t *itf, size_t N, uint64_t m, uint64_t invN, const uint64_t *prec_itf, uint64_t prec_invN) {
    size_t t = 1;
    size_t log1 = 1;
    size_t h, j, i, j1, j2;
    uint64_t u, v, vv, s, s_prec;
    for (h = (N>>1); h >= 1; h >>= 1) {
        for (i = 0; i < h; i++) {
            j1 = i << log1;
            j2 = j1 + t;
            s = itf[h + i]; 
            s_prec = prec_itf[h + i];
            for (j = j1; j < j2; j++) {
                u = ax[j];
                v = ax[j + t];

                vv = u + v;
                ax[j] = vv < m ? vv : vv - m;

                v = (u > v) ? u - v : m + u - v;

                ModMulShoupEq(v, s, m, s_prec);

                ax[j+t] = v;
            }
        }
        t <<= 1;
        log1++;
    }

    for (size_t j = 0; j < N; j++) {
        ModMulShoupEq(ax[j], invN, m, prec_invN);
    }
}

}