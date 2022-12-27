#include <iostream>

#include "src/include/math.h"

using namespace std;

int main(int argc, char const *argv[])
{
    uint64_t m = 17;
    uint64_t logm = 5;
    uint64_t a = 9;
    uint64_t b = 11;

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "m = " << m << endl;

    uint64_t z = BarrettPrecompute(m, logm);
    cout << "z = " << z << endl;

    cout << "a + b mod m = " << ModAdd(a, b, m) << endl;
    cout << "a * b mod m = " << ModMultBarrett(a, b, m, z, logm) << endl;
    cout << "(check) a * b mod m = " << (a * b) % m << endl;


    uint64_t s = ShoupPrecompute(b, m);
    cout << "a * b mod m = " << ModMulShoup(a, b, m, s) << endl;
    
    return 0;
}

/*
пример больших чисел для большого умножения...
    uint64_t m = 4294967311;
    uint64_t logm = 33;
    uint64_t a = 4294967311 - 8239042;
    uint64_t b = 4294967311 - 23984778;
*/
