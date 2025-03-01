#include <iostream>
#include <string>
#include <iomanip>

#include "vectbigint.h"
#include "examples_utils.h"

using namespace std;
using namespace NTL;
using namespace polycalc;


int main(int argc, char const *argv[])
{
    uint64_t q = (1ULL << 30) - (1ULL << 18) + 1;
    size_t logN = 10;
    size_t N = 1 << logN;

    size_t s = 8;
    size_t bsize = s*N / 2;
    size_t r_bits = bsize;
    uint64_t mu;


    if (argc > 1) {
        int r = stoi(argv[1]); 
        bsize = size_t(1) << static_cast<size_t>(ceil(log2(double(r))));
    }

    if (argc > 2) {
        s = stoi(argv[3]); 
    }
    

    NTTInstance ntt = NTTManager::GetNTTPtr(N, q);

    uint64_t mx[N];
    uint64_t rsx[N];

    uint64_t ax[N];
    uint64_t bx[N];
    uint64_t cx[N];

    ZZ a, b, c, c_ref;
    ZZ a_bar, b_bar, c_bar;
    ZZ m;
    m = GenPrime_ZZ(bsize);

    REDC_precompute(mx, rsx, mu, m, r_bits, s, N);


    cout << "Parameters: " << endl;
    cout << setw(8) << "q = "     << q     << endl; 
    cout << setw(8) << "logN = "  << logN  << endl; 
    cout << setw(8) << "N = "     << N     << endl; 
    cout << setw(8) << "bsize = " << bsize << endl; 
    cout << setw(8) << "s = "     << s     << endl; 
    cout << setw(8) << "mu = "    << mu    << endl; 
    // PrintAll(mx, N, s, "m");

    a = RandomBits_ZZ(bsize) % m;
    b = RandomBits_ZZ(bsize) % m;


    // a_bar = MontIn(a, r_bits, m);
    // b_bar = MontIn(b, r_bits, m);
    // cout << setw(8) << "a_bar = "     << a_bar  << endl; 
    // cout << setw(8) << "b_bar = "     << b_bar  << endl; 


    TransformZZtoCPoly(ax, a, N, s);
    TransformZZtoCPoly(bx, b, N, s);

    
    c_ref = (a * b) % m;

    cout << setw(8) << "m = "     << m     << endl; 
    cout << setw(8) << "a = "     << a     << endl; 
    cout << setw(8) << "b = "     << b     << endl; 
    cout << setw(8) << "c_ref = " << c_ref << endl; 

    REDC_In(ax, rsx, mx, mu, r_bits, s, N, ntt);
    REDC_In(bx, rsx, mx, mu, r_bits, s, N, ntt);



    QMul(cx, ax, bx, ntt);
    CarryPropagation(cx, N, s);
    REDC(cx, mx, mu, r_bits, s, N);

    REDC_Out(cx, mx, mu, r_bits, s, N);

    ReconstructZZfromCPoly(c, cx, N, s);

    //c = MontOut(c_bar, r_bits, m);

    cout << setw(8) << "c = " << c << endl; 
    cout << setw(8) << "diff = " << c - c_ref << endl; 

    return 0;
}


