#include <iostream>

#include "math.h"
#include "polynomial.h"

using namespace std;

int main(int argc, char const *argv[])
{
    uint64_t m = 17;
    uint64_t logm = MSB(m);
    uint64_t a = 9;
    uint64_t b = 11;


//пример больших чисел для большого умножения...
    // uint64_t m = 4294967311;
    // uint64_t logm = MSB(m);
    // uint64_t a = 4294967311 - 8239042;
    // uint64_t b = 4294967311 - 23984778;


    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "logm = " << logm << endl;
    cout << "m = " << m << endl;

    uint64_t z = BarrettPrecompute(m, logm);
    cout << "z = " << z << endl;

    cout << "a + b mod m = " << ModAdd(a, b, m) << endl;
    cout << "a * b mod m = " << ModMultBarrett(a, b, m, z, logm) << endl;
    cout << "(check) a * b mod m = " << (a * b) % m << endl;


    uint64_t s = ShoupPrecompute(b, m);
    cout << "a * b mod m = " << ModMulShoup(a, b, m, s) << endl;


    Polynomial poly(8, 11);

    uint64_t sum = 0;
    for (size_t i = 0; i < 8; i++) {
        poly[i] = i+1;
        sum += i+1;
    }

    cout << "A(x) = " << poly << endl;

    cout << poly(1) << endl;
    cout << sum % poly.GetModulus() << endl;

    cout << poly(0) << endl;
    cout << poly(3) << endl;

    Polynomial poly1(8, 11);
    for (size_t i = 0; i < 8; i++) {
        poly1[i] = (i*i) % 11;
    }
    cout << "B(x) = "  << poly1 << endl;

    Polynomial poly2;
    poly2 = poly + poly1;

    cout << "A(x) + B(x) = " << poly2 << endl;
    
    return 0;
}


