#include <iostream>

#include "ntmath.h"

using namespace std;

int main(int argc, char const *argv[]) {
    size_t logN = 11;
    size_t logm = 58;
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


    uint64_t a = (uint64_t(1L) << logm) - 1;

    cout << "a = " << a << endl;
    uint64_t b = RhoPollard(a);
    cout << "b = " << b << endl;

    vector<uint64_t> fact = Factorize(a);

    cout << "factorize(a) = ";
    uint64_t check_a = 1;
    for (size_t i = 0; i < fact.size(); i++) {
        cout << fact[i] << " ";
        check_a *= fact[i];
    }
    cout << endl;
    cout << "check_a = " << check_a << endl;

    vector<uint64_t> factprime = Factorize(m);

    cout << "m = " << m << endl;
    cout << "factorize(m) = ";
    for (size_t i = 0; i < factprime.size(); i++) {
        cout << factprime[i] << " ";
    }
    cout << endl;
    
    uint64_t g = FindPrimitive(m);
    uint64_t gM = FindGenerator(m, M);
    cout << "m = " << m << endl;
    cout << "g = " << g << endl;
    cout << "gM = " << gM << endl;

    cout << "g^phi(m) = " << ModExp(g, m - 1, m) << endl;

    for (size_t i = 1; i <= M; i<<=1) {
        cout << "gM^(" << i << ") = " << ModExp(gM, i, m) << endl;
    }

    return 0;
}