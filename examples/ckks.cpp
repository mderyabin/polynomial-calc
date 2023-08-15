#include <iostream>
#include "encode.h"

#include "examples_utils.h"

using namespace std;
using namespace polycalc;

int main(int argc, char const *argv[]) {
    // only independent parameters, specified by user
    size_t logN = 8;
    size_t logq = 20;

    size_t slots = 2;
    size_t scale = 1 << 6;

    // read parameters from command line (logN first, logm second)
    read_command_line(logN, logq, argc, argv);

    size_t N = 1 << logN;
    size_t M = N << 1;

    cout << "slots = " << slots << endl;
    cout << "scale = " << scale << endl;
    cout << "logq = " << logq << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;

    uint64_t q = FindFirstPrimeUp(logq, M);
    cout << "q = " << q << endl;


    vector<complex<long double>> a(slots);
    for (size_t i = 0; i < slots; i++) {
        a[i] = complex<long double>(generaternd(-5.0,5.0),  generaternd(-5.0,5.0));
    }

    cout << "input a: ";
    for (size_t i = 0; i < a.size(); i++)
        cout << a[i] << " ";
    cout << endl;

    Polynomial ra = encode(a, N, q, scale);

    cout << "encoded polynomial for a: ";
    cout << ra << endl;

    vector<complex<long double>> aa = decode(ra, scale);

    cout << "a check: ";
    for (size_t i = 0; i < aa.size(); i++)
        cout << aa[i] << " ";
    cout << endl;

    return 0;
}


