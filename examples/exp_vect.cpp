#include <iostream>
#include <iomanip>

#include "vectbigint.h"
#include "examples_utils.h"

using namespace std;
using namespace NTL;
using namespace polycalc;


ZZ FastModExp(ZZ a, ZZ e, ZZ m);

int main(int argc, char const *argv[])
{
    uint64_t q = (1ULL << 30) - (1ULL << 18) + 1;
    size_t logN = 8;
    size_t N = 1 << logN;

    size_t s = 8;
    size_t bsize = s*N / 2;
    size_t r_bits = bsize;
    uint64_t mu;

    NTTInstance ntt = NTTManager::GetNTTPtr(N, q);

    uint64_t mx[N];
    uint64_t rsx[N];

    uint64_t ax[N];
    uint64_t bx[N];
    uint64_t cx[N];

    ZZ a, b, c, c_ref;
    ZZ a_bar, b_bar, c_bar;
    ZZ m;
    m = GenPrime_ZZ(bsize-2);

    REDC_precompute(mx, rsx, mu, m, r_bits, s, N);

    cout << "Parameters: " << endl;
    cout << setw(8) << "q = "     << q     << endl; 
    cout << setw(8) << "logN = "  << logN  << endl; 
    cout << setw(8) << "N = "     << N     << endl; 
    cout << setw(8) << "bsize = " << bsize << endl; 
    cout << setw(8) << "s = "     << s     << endl; 
    cout << setw(8) << "mu = "    << mu    << endl; 

    a = RandomBits_ZZ(bsize) % m;
    b = RandomBits_ZZ(bsize) % m;

    TransformZZtoCPoly(ax, a, N, s);
    TransformZZtoCPoly(bx, b, N, s);

    cout << setw(8) << "m = "     << m     << endl; 
    cout << setw(8) << "a = "     << a     << endl; 
    cout << setw(8) << "b = "     << b     << endl;


    cout << " -- exp mod test -- " << endl;
    FastModExp(cx, ax, bx, mx, N, s);
    //cout << " cx = ax ^ bx mod mx = "; PrintArray(cx, N);
    cout << "         a ^ b mod m = " << FastModExp(a, b, m) << endl;
    ReconstructZZfromCPoly(c, cx, N, s);
    cout << "     c = a ^ b mod m = " << c << endl;
    cout << "                diff = " << c-FastModExp(a, b, m) << endl;

    FastModExpREDC(cx, ax, bx, mx, rsx, mu, r_bits, N, s, ntt);
    ReconstructZZfromCPoly(c, cx, N, s);
    cout << "REDC c = a ^ b mod m = " << c << endl;
    cout << "                diff = " << c-FastModExp(a, b, m) << endl;

    return 0;
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