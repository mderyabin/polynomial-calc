#include <iostream>

#include "polynomial.h"

using namespace std;

int main(int argc, char const *argv[])
{
    // uint64_t m = 4294967311;
    uint64_t m = 17;
	uint64_t N = 8;
	uint64_t logm = MSB(m);

    Polynomial poly(N, m);

    uint64_t sum = 0;
    for (size_t i = 0; i < 8; i++) {
        poly[i] = i+1;
        sum += i+1;
    }

    cout << "A(x) = " << poly << endl;

    cout << "A(0) = " << poly(1) << endl;
    cout << "check " << sum % poly.GetModulus() << endl;

    cout << "A(0) = " << poly(0) << endl;
    cout << "A(3) = " << poly(3) << endl;

    Polynomial poly1(N, m);
    for (size_t i = 0; i < 8; i++) {
        poly1[i] = (i*i) % 11;
    }
    cout << "B(x) = "  << poly1 << endl;

    Polynomial poly2;
    poly2 = poly + poly1;

    cout << "A(x) + B(x) = " << poly2 << endl;

    Polynomial poly3;
    poly3 = poly * poly1;

    cout << "A(x) * B(x) mod (X^"<< poly.GetN() <<" + 1) = " << poly3 << endl;

	Polynomial monomial(N, m, true);
    monomial[2] = 1;
	
	cout << "D(X) = " << monomial << endl;
	cout << "A(X) * D(X) = " << poly * monomial << endl;
	
	
    return 0;
}


