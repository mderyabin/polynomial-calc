#include <iostream>
#include "encode.h"
#include "ntmath.h"

#include "examples_utils.h"

using namespace std;
using namespace polycalc;

int main(int argc, char const *argv[]) {
    // only independent parameters, specified by user
    size_t logN = 3;
    size_t logq = 7;

    // read parameters from command line (logN first, logq second)
    read_command_line(logN, logq, argc, argv);

    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeUp(logq, M);

    cout << "logq = " << logq << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;
    cout << "q = " << q << endl;

    uint64_t m_raw[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Polynomial m(m_raw, N, q);
    cout << "m(x) = " << m << endl;

    Polynomial a(N, q);
    a.GenerateUniform();
    cout << "a(x) = " << a << endl;

    Polynomial e(N, q);
    e.GenerateDiscreteGauss();
    cout << "e(x) = " << e << endl;

    uint64_t delta = 1 << 4;
    cout << "delta = " << delta << endl;

    Polynomial s(N, q);
    s.GenerateBinary();
    cout << "s(x) = " << s << endl;

    // (X^8 + 1)
    // X^8 = -1

    Polynomial dm(N, q);
    for (size_t i = 0; i < N; i++) {
        dm[i] = (m[i] * delta) % q;
    }
    cout << "dm(x) = " << dm << endl;

    Polynomial b(N, q);
    b = a*s + dm + e;

    cout << "b(x) = " << b << endl;

    
    Polynomial dec_dm(N, q);

    dec_dm = b - a*s;
    cout << "dec_dm(x) = " << dec_dm << endl;

    Polynomial dec_m(N, q);

    for (size_t i = 0; i < N; i++) {
        dec_m[i] = static_cast<uint64_t>(round( static_cast<double>(dec_dm[i]) / (1.0*delta) )); 
    }
    
    cout << "dec_m(x) = " << dec_m << endl;

    return 0;
}