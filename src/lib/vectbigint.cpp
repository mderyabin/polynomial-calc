#include "vectbigint.h"

using namespace std;
using namespace NTL;

namespace polycalc {

void LeftShift(uint64_t *cx, const uint64_t *ax, size_t shift, size_t n, size_t s) {
    size_t bs = shift % s;
    size_t bc = shift / s;

    // cout << "block_count = " << bc << endl;
    // cout << "bit_shift = " << bs << endl;
    
    size_t mask = (1 << s) - 1;
    uint64_t clo, chi;

    cx[bc] = (ax[0] << bs) & mask;
    for (size_t i = 1; i < n-bc; i++) {
        clo = (ax[i] << bs) & mask;
        chi = ax[i-1] >> (s - bs);
        cx[bc + i] = clo + chi;
    }
}

void LeftShift(uint64_t *cx, size_t shift, size_t n, size_t s) {

    uint64_t *ax = new uint64_t[n];
    // uint64_t *ax = globalMemPool.tmp2;

    size_t bs = shift % s;
    size_t bc = shift / s;

    copy(cx, cx+n, ax);

    // cout << "block_count = " << bc << endl;
    // cout << "bit_shift = " << bs << endl;
    
    size_t mask = (1 << s) - 1;
    uint64_t clo, chi;

    cx[bc] = (ax[0] << bs) & mask;
    for (size_t i = 1; i < n-bc; i++) {
        clo = (ax[i] << bs) & mask;
        chi = ax[i-1] >> (s - bs);
        cx[bc + i] = clo + chi;
    }

    delete [] ax;
}

bool IsBiggerOrEqual(const uint64_t *ax, const uint64_t *bx, size_t n) {
    size_t i = n-1;
    while (ax[i] == bx[i] and (i < n)) i--;
    if (ax[i] < bx[i]) return false;
    else return true;
}

bool IsBigger(const uint64_t *ax, const uint64_t *bx, size_t n) {
    size_t i = n-1;
    while (ax[i] == bx[i] and (i < n)) i--;
    if (ax[i] > bx[i]) return true;
    else return false;
}


uint64_t CarryPropagation(uint64_t *cx, const uint64_t *ax, size_t n, size_t s) {
    uint64_t carry = 0, cxi;
    size_t mask = (1 << s) - 1;
    for (size_t i = 0; i < n; i++) {
        cxi = ax[i] + carry;
        cx[i] = cxi & mask;
        carry = cxi >> s;
    }
    return carry;
}

uint64_t CarryPropagation(uint64_t *cx, size_t n, size_t s) {
    uint64_t carry = 0, cxi;
    size_t mask = (1 << s) - 1;
    for (size_t i = 0; i < n; i++) {
        cxi = cx[i] + carry;
        cx[i] = cxi & mask;
        carry = cxi >> s;
    }
    return carry;
}

void AddLazy(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t n) {
    for (size_t i = 0; i < n; i++) {
        cx[i] = ax[i] + bx[i];
    }
}

void AddLazy(uint64_t *cx, const uint64_t *ax, size_t n) {
    for (size_t i = 0; i < n; i++) {
        cx[i] = cx[i] + ax[i];
    }
}

uint64_t AddFull(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t n, size_t s) {
    AddLazy(cx, ax, bx, n);
    return CarryPropagation(cx, n, s);
}

uint64_t AddFull(uint64_t *cx, const uint64_t *ax, size_t n, size_t s) {
    AddLazy(cx, ax, n);
    return CarryPropagation(cx, n, s);
}

// a > b!!
uint64_t SubLazy(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t n, size_t s) {
    uint64_t neg_carry = 0;
    uint64_t axi, bxi;
    uint64_t base = 1 << s;
    for (size_t i = 0; i < n; i++) {
        axi = ax[i];
        bxi = bx[i]+neg_carry;
        if (axi < bxi) {
            axi += base;
            neg_carry = 1;
        } else {
            neg_carry = 0;
        }
        cx[i] = axi - bxi;
    }
    return neg_carry;
}

uint64_t SubLazy(uint64_t *cx, const uint64_t *bx, size_t n, size_t s) {
    uint64_t neg_carry = 0;
    uint64_t cxi, bxi;
    uint64_t base = 1 << s;
    for (size_t i = 0; i < n; i++) {
        cxi = cx[i];
        bxi = bx[i] + neg_carry;
        if (cxi < bxi) {
            cxi += base;
            neg_carry = 1;
        } else {
            neg_carry = 0;
        }
        cx[i] = cxi - bxi;
    }
    return neg_carry;
}

void AddMod(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s) {
    AddFull(cx, ax, bx, n, s);
    if (IsBiggerOrEqual(cx, mx, n)) {
        SubLazy(cx, mx, n, s);
    }
}

void AddMod(uint64_t *cx, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s) {
    AddFull(cx, bx, n, s);
    if (IsBiggerOrEqual(cx, mx, n)) {
        SubLazy(cx, mx, n, s);
    }
}

void AddModLazy(uint64_t *cx, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s) {
    AddLazy(cx, bx, n);
    if (IsBiggerOrEqual(cx, mx, n)) {
        SubLazy(cx, mx, n, s);
    }
}

void MulMod(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s) {
    uint64_t *shifted = new uint64_t[n];
    // uint64_t *shifted = globalMemPool.shifted;

    copy(ax, ax+n, shifted);

    if (bx[0] % 2 == 1) {
        copy(ax, ax+n, cx);
    } else {
        for (size_t i = 0; i < n; i++) {
            cx[i] = 0;
        }
    }
   
    for (size_t i = 1; i < n*s; i++) {
        LeftShift(shifted, 1, n, s);

        while (IsBiggerOrEqual(shifted, mx, n)) {
            SubLazy(shifted, mx, n, s);
        }
        
        size_t j = i/s;
        size_t k = i%s;

        if ((bx[j] >> k)%2 == 1) AddMod(cx, shifted, mx, n, s);

    }

    // CarryPropagation(cx, n, s);


    delete [] shifted;
}

void MulMod(uint64_t *cx, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s) {
    uint64_t *shifted = new uint64_t[n];
    // uint64_t *shifted = globalMemPool.shifted;

    copy(cx, cx+n, shifted);

    if (bx[0] % 2 == 0) {
        for (size_t i = 0; i < n; i++) {
            cx[i] = 0;
        }
    }
   
    for (size_t i = 1; i < n*s; i++) {
        LeftShift(shifted, 1, n, s);

        while (IsBiggerOrEqual(shifted, mx, n)) {
            SubLazy(shifted, mx, n, s);
        }
        
        size_t j = i/s;
        size_t k = i%s;

        if ((bx[j] >> k)%2 == 1) AddMod(cx, shifted, mx, n, s);

    }

    // CarryPropagation(cx, n, s);

    delete [] shifted;
}

void SqrMod(uint64_t *cx, const uint64_t *ax, const uint64_t *mx, size_t n, size_t s) {
    // copy(ax, ax+n, cx);
    MulMod(cx, ax, ax, mx, n, s);
}

void SqrMod(uint64_t *cx, const uint64_t *mx, size_t n, size_t s) {
    uint64_t *shifted = new uint64_t[n];
    uint64_t *bx = new uint64_t[n];

    // uint64_t *shifted = globalMemPool.shifted;
    // uint64_t *bx = globalMemPool.tmp1;

    copy(cx, cx+n, shifted);
    copy(cx, cx+n, bx);  

    if (bx[0] % 2 == 0) {
        for (size_t i = 0; i < n; i++) {
            cx[i] = 0;
        }
    } 
   
    for (size_t i = 1; i < n*s; i++) {
        LeftShift(shifted, 1, n, s);

        while (IsBiggerOrEqual(shifted, mx, n)) {
            SubLazy(shifted, mx, n, s);
        }
        
        size_t j = i/s;
        size_t k = i%s;

        if ((bx[j] >> k)%2 == 1) AddMod(cx, shifted, mx, n, s);

    }

    // CarryPropagation(cx, n, s);

    delete [] bx;
    delete [] shifted;
}

void FastModExp(uint64_t *cx, const uint64_t *ax, const uint64_t *ex, const uint64_t *mx, size_t n, size_t s) { 

    uint64_t *base = new uint64_t[n];
    // uint64_t *base = globalMemPool.base;

    copy(ax, ax+n, base);

    bool flag = false;

    size_t emsb = 0;
    for (size_t i = 0; i < n; i++) {
        if (ex[i] != 0) emsb = i;
    }
    emsb++;

    for (size_t i = 0; i < emsb*s; i++) {
        size_t j = i/s;
        size_t k = i%s;


        if ((ex[j] >> k)%2 == 1) {
            if (flag) { 
                MulMod(cx, base, mx, n, s);
            } else {
                copy(base, base+n, cx);
                flag = true;
            }
        }

        if (! (j >= emsb and (ex[j] >> k) == 0) ) {
            SqrMod(base, mx, n, s);
        }
    }

    delete [] base;
}

uint64_t MulLazy(uint64_t *cx, const uint64_t *ax, uint64_t b, size_t n) {
    for (size_t i = 0; i < n; i++) {
        cx[i] = ax[i] * b;
    }
}

uint64_t MulLazy(uint64_t *cx, uint64_t b, size_t n) {
    for (size_t i = 0; i < n; i++) {
        cx[i] *= b;
    }
}

void QMul(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, NTTInstance ntt) {
    size_t N = ntt->GetN();
    uint64_t m = ntt->GetModulus();
    size_t logm = ntt->GetLogModulus();
    uint64_t prec = ntt->GetPrecBarrett();

    uint64_t ax_ntt[N];
    uint64_t bx_ntt[N];
    uint64_t cx_ntt[N];

    ntt->ComputeForward(ax_ntt, ax);
    ntt->ComputeForward(bx_ntt, bx);

    ModHadamardMul(cx_ntt, ax_ntt, bx_ntt, N, m, prec, logm);

    ntt->ComputeInverse(cx, cx_ntt);
}

void QMul(uint64_t *cx, const uint64_t *bx, NTTInstance ntt) {
    size_t N = ntt->GetN();
    uint64_t m = ntt->GetModulus();
    size_t logm = ntt->GetLogModulus();
    uint64_t prec = ntt->GetPrecBarrett();

    uint64_t ax_ntt[N];
    uint64_t bx_ntt[N];
    uint64_t cx_ntt[N];

    ntt->ComputeForward(ax_ntt, cx);
    ntt->ComputeForward(bx_ntt, bx);

    ModHadamardMul(cx_ntt, ax_ntt, bx_ntt, N, m, prec, logm);

    ntt->ComputeInverse(cx, cx_ntt);
}

void REDC_In(uint64_t *cx, const uint64_t *ax, const uint64_t *rsx, const uint64_t *mx, uint64_t mu, size_t r_bits, size_t s, size_t N, NTTInstance ntt) {
    QMul(cx, ax, rsx, ntt);
    REDC(cx, mx, mu, r_bits, s, N);
}

void REDC_In(uint64_t *cx, const uint64_t *rsx, const uint64_t *mx, uint64_t mu, size_t r_bits, size_t s, size_t N, NTTInstance ntt) {
    QMul(cx, rsx, ntt);
    REDC(cx, mx, mu, r_bits, s, N);
}

void REDC_Out(uint64_t *cx, const uint64_t *mx, uint64_t mu, size_t r_bits, size_t s, size_t N) {
    REDC(cx, mx, mu, r_bits, s, N);
}

void REDC_Out(uint64_t *cx, const uint64_t *ax, const uint64_t *mx, uint64_t mu, size_t r_bits, size_t s, size_t N) {
    copy(ax, ax+N, cx);
    REDC(cx, mx, mu, r_bits, s, N);
}

void REDC(uint64_t *cx, const uint64_t *mx, uint64_t mu, size_t r_bits, size_t s, size_t N) {
    size_t r_blocks = r_bits / s;
    uint64_t b_mask = (1ULL << s) - 1;

    uint64_t tx[N];
    uint64_t tmp[N];

    copy(cx, cx+N, tx);

    uint64_t carry = 0;
    uint64_t x, qi;

    for (size_t i = 0; i < r_blocks; i++) {
        // carry = 0;
        qi = (tx[i] * mu) & b_mask;

        MulLazy(tmp, mx, qi, N - r_blocks);
        AddLazy(tx + i, tmp, N - r_blocks);
        CarryPropagation(tx + i, N-i, s);

        // for (size_t j = 0; j < N - r_blocks; j++) {
        //     x = tx[i+j] + qi*mx[j] + carry;
        //     tx[i+j] = x & b_mask;
        //     carry = (x >> s);
        // }

        // for (size_t j = N - r_blocks; j < N - i; j++) {
        //     x = tx[i+j] + carry;
        //     tx[i+j] = x & b_mask;
        //     carry = (x >> s);
        // }
    }

    fill(cx, cx+N, 0);
    for (size_t i = 0; i < N - r_blocks; i++) {
        cx[i] = tx[r_blocks + i];
    } 
}

void FastModExpREDC(uint64_t *cx, const uint64_t *ax, const uint64_t *ex, const uint64_t *mx, const uint64_t *rsx, uint64_t mu, size_t r_bits, size_t n, size_t s, NTTInstance ntt) { 

    uint64_t *base = new uint64_t[n];
    // uint64_t *base = globalMemPool.base;

    copy(ax, ax+n, base);
    REDC_In(base, rsx, mx, mu, r_bits, s, n, ntt);

    bool flag = false;

    size_t emsb = 0;
    for (size_t i = 0; i < n; i++) {
        if (ex[i] != 0) emsb = i;
    }
    emsb++;

    for (size_t i = 0; i < emsb*s; i++) {
        size_t j = i/s;
        size_t k = i%s;


        if ((ex[j] >> k)%2 == 1) {
            if (flag) { 
                QMul(cx, base, ntt);
                CarryPropagation(cx, n, s);
                REDC(cx, mx, mu, r_bits, s, n);
            } else {
                copy(base, base+n, cx);
                flag = true;
            }
        }

        if (! (j >= emsb and (ex[j] >> k) == 0) ) {
            QMul(base, base, base, ntt);
            CarryPropagation(base, n, s);
            REDC(base, mx, mu, r_bits, s, n);
        }
    }

    REDC_Out(cx, mx, mu, r_bits, s, n);

    delete [] base;
}

void TransformZZtoCPoly(uint64_t *ax, const NTL::ZZ a, size_t n, size_t s) {    
    NTL::ZZ a_copy = a;
    for (size_t i = 0; i < n; i++) {
        ax[i] = a_copy % (1<<s);
        a_copy /= (1<<s);
    }
}

void REDC_precompute(uint64_t *mx, uint64_t *rsx, uint64_t &mu, ZZ m, size_t r_bits, size_t s, size_t N) {
    ZZ r, r_sqr;
    r = (ZZ(1) << r_bits) % m;
    r_sqr = (r << r_bits) % m;
    // InvMod(r_inv, r, m);

    TransformZZtoCPoly(rsx, r_sqr, N, s);
    TransformZZtoCPoly(mx, m, N, s);


    uint64_t b = 1ULL << s;
    uint64_t m_mod_b = m % b;
    uint64_t m_inv_b;

    m_inv_b = InvMod(m_mod_b, b);

    mu = b - m_inv_b;

}

// void MemoryPool::Init(size_t n) {
//     base = new uint64_t[n];
//     shifted = new uint64_t[n];
//     tmp1 = new uint64_t[n];
//     tmp2 = new uint64_t[n];
//     tmp3 = new uint64_t[n];
// }

// MemoryPool::~MemoryPool() {
//     delete [] tmp3;
//     delete [] tmp2;
//     delete [] tmp1;
//     delete [] shifted;
//     delete [] base;
// }

}