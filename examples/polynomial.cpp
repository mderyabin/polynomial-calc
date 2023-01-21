#include <iostream>

#include "polynomial.h"

using namespace std;

int main(int argc, char const *argv[])
{
    size_t logN = 3;
    size_t logm = 5;
    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t m = FindFirstPrimeUp(logm, M);

    cout << "logm = " << logm << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;
    cout << "m = " << m << endl;

    Polynomial A(N, m);
    A.GenerateUniform();

    Polynomial B(N, m);
    B.GenerateUniform();

    cout << "A(x) = " << A << endl;
    cout << "B(x) = " << B << endl;

    cout << "naive A(x) * B(x) = " << A*B << endl;

    A.SetFormatEval();
    B.SetFormatEval();

    cout << "NTT(A(x)) = " << A << endl;
    cout << "NTT(B(x)) = " << B << endl;

    Polynomial C = A * B;

    cout << "NTT(C) = NTT(A(x)) * NTT(B(x)) = " << C << endl;

    C.SetFormatCoef();

    cout << "C = INTT(NTT(C)) = " << C << endl;
	
    return 0;
}


