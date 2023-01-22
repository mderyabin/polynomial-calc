#include "utils.h"

#include "encode.h"

#include <random>
#include <vector>
#include <complex>

using namespace std;
using namespace polycalc;

long double generaternd(long double a, long double b);

int main(int argc, char const *argv[]) {

    // only independent parameters, specified by user
    size_t logN = 4;
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

    vector<complex<long double>> b(slots);
    for (size_t i = 0; i < slots; i++) {
        b[i] = complex<long double>(generaternd(-5.0,5.0),  generaternd(-5.0,5.0));
    }

    cout << "input b: ";
    for (size_t i = 0; i < b.size(); i++)
        cout << b[i] << " ";
    cout << endl;

    Polynomial rb = encode(b, N, q, scale);

    cout << "encoded polynomial for a: ";
    cout << ra << endl;
    cout << "encoded polynomial for b: ";
    cout << rb << endl;

    ra.SetFormatEval();
    rb.SetFormatEval();

    Polynomial rc = ra * rb;

    rc.SetFormatCoef();

    cout << "polynomial for c = a*b: ";
    cout << rc << endl;

    vector<complex<long double>> c = decode(rc, scale*scale);

    cout << "c = a*b after decode: ";
    for (size_t i = 0; i < slots; i++)
        cout << c[i] << " ";
    cout << endl;

    cout << "a*b check: ";
    for (size_t i = 0; i < a.size(); i++)
        cout << a[i]*b[i] << " ";
    cout << endl;

    return 0;
}

long double generaternd(long double a, long double b) {
    random_device rd;
    uniform_real_distribution<long double> dist(a, b);
    default_random_engine en(rd());
    return dist(en);
}