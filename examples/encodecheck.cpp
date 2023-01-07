#include <iostream>

#include "polynomial.h"
#include "encode.h"

#include <vector>
#include <complex>

using namespace std;

int main(int argc, char const *argv[]) {
    uint64_t m = 2048;
    size_t N = 4;

    vector<complex<long double>> a(2);
    a[0] = complex<long double>(3.0,  4.0);
    a[1] = complex<long double>(2.0, -1.0);

    cout << "input: ";
    for (size_t i = 0; i < a.size(); i++)
        cout << a[i] << " ";
    cout << endl;
    Polynomial r = encode(a, N, m, 128);

    cout << r << endl;

    vector<complex<long double>> b = decode(r, 128);

    cout << "input after decode: ";
    for (size_t i = 0; i < b.size(); i++)
        cout << b[i] << " ";
    cout << endl;

    return 0;
}
