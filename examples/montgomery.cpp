#include <iostream>
#include <string>
#include <iomanip>

#include "vectbigint.h"
#include "examples_utils.h"

using namespace std;
using namespace NTL;
using namespace polycalc;

ZZ MontIn(ZZ a, size_t r_bits, ZZ m);
ZZ MontOut(ZZ a_bar, size_t r_bits, ZZ m);

void REDC(uint64_t *cx, const uint64_t *mx, size_t s, size_t N);

int main(int argc, char const *argv[])
{
    uint64_t q = (1ULL << 30) - (1ULL << 18) + 1;
    size_t logN = 8;
    size_t N = 1 << logN;

    size_t s = 8;
    size_t bsize = s*N / 2;
    size_t r_bits = bsize;


    if (argc > 1) {
        int r = stoi(argv[1]); 
        bsize = size_t(1) << static_cast<size_t>(ceil(log2(double(r))));
    }

    if (argc > 2) {
        s = stoi(argv[3]); 
    }

    cout << "Parameters: " << endl;
    cout << setw(8) << "q = "     << q     << endl; 
    cout << setw(8) << "logN = "  << logN  << endl; 
    cout << setw(8) << "N = "     << N     << endl; 
    cout << setw(8) << "bsize = " << bsize << endl; 
    cout << setw(8) << "s = "     << s     << endl; 
    

    NTTInstance ntt = NTTManager::GetNTTPtr(N, q);

    uint64_t mx[N];

    uint64_t ax[N];
    uint64_t bx[N];
    uint64_t cx[N];

    ZZ a, b, c, c_ref;
    ZZ a_bar, b_bar, c_bar;
    ZZ m;
    m = GenPrime_ZZ(bsize);

    a = RandomBits_ZZ(bsize) % m;
    b = RandomBits_ZZ(bsize) % m;

    TransformZZtoCPoly(mx, m, N, s);

    a_bar = MontIn(a, r_bits, m);
    b_bar = MontIn(b, r_bits, m);

    TransformZZtoCPoly(ax, a_bar, N, s);
    TransformZZtoCPoly(bx, b_bar, N, s);

    c_ref = (a * b) % m;

    cout << setw(8) << "m = "     << m     << endl; 
    cout << setw(8) << "a = "     << a     << endl; 
    cout << setw(8) << "b = "     << b     << endl; 
    cout << setw(8) << "c_ref = " << c_ref << endl; 

    QMul(cx, ax, bx, ntt);
    CarryPropagation(cx, N, s);
    
    REDC(cx, mx, s, N);

    ReconstructZZfromCPoly(c_bar, cx, N, s);

    c = MontOut(c_bar, r_bits, m);

    cout << setw(8) << "c = " << c << endl; 
    cout << setw(8) << "diff = " << c - c_ref << endl; 

    return 0;
}

ZZ MontIn(ZZ a, size_t r_bits, ZZ m) {
    ZZ a_bar = (a << r_bits) % m;
    return a_bar;
}

ZZ MontOut(ZZ a_bar, size_t r_bits, ZZ m) {
    ZZ r, r_inv;
    r = (ZZ(1) << r_bits) % m;
    InvMod(r_inv, r, m);
    ZZ a = (a_bar * r_inv) % m;
    return a;
}

void REDC(uint64_t *cx, const uint64_t *mx, size_t s, size_t N) {
    ZZ m, r;

    size_t r_bits;
    size_t r_blocks;

    uint64_t b;
    uint64_t m_prime;
    uint64_t b_mask = (1ULL << s) - 1;

    uint64_t tx[N];
    uint64_t tmp[N];
    copy(cx, cx+N, tx);

    r_blocks = N/2;
    r_bits = s * r_blocks;
    b = 1ULL << s;

    ReconstructZZfromCPoly(m, mx, r_blocks, s);
    

    

    uint64_t m_mod_b = m % b;
    uint64_t m_inv_b;

    m_inv_b = InvMod(m_mod_b, b);

    m_prime = b - m_inv_b;


    uint64_t carry = 0;
    uint64_t x, qi;

    for (size_t i = 0; i < r_blocks; i++) {
        // carry = 0;
        qi = (tx[i] * m_prime) & b_mask;

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
