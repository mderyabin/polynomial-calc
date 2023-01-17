#include <iostream>

#include "ntmath.h"

using namespace std;

int main(int argc, char const *argv[]) {
    size_t logN = 11;
    size_t logm = 40;
    size_t N = 1 << logN;
    size_t M = N << 1;

    cout << "logm = " << logm << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;

    uint64_t m = FindFirstPrimeDown(logm, M);
    cout << "m[1] = " << m << endl;

    m = FindPrevPrime(m, M);
    cout << "m[2] = " << m << endl;

    return 0;
}