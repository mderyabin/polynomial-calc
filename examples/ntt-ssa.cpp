// https://dl.acm.org/doi/pdf/10.1145/1277548.1277572

#include <iostream>
#include "encode.h"
#include "ntmath.h"
#include "polynomial.h"

#include <NTL/ZZ.h>

#include "examples_utils.h"

using namespace std;
using namespace polycalc;
using namespace NTL;

const size_t osize_bit = 1024;

const size_t N = 2*osize_bit;       // size of numbers (product of a and b size) 
const size_t w = 64;                // word size

const size_t M = 16;                 // coefficients size
const size_t k = N/M;               // number of coefficients

const size_t log_q = 58;            // size of q

void TransformZZtoCPoly(Polynomial &ax, const ZZ a);
void ReconstructZZfromCPoly(ZZ &a, const Polynomial ax);

int main() {

    cout << "osize_bit = " << osize_bit << endl;
    cout << "N = " << N << endl;
    cout << "w = " << w << endl;
    cout << "M = " << M << endl;
    cout << "k = " << k << endl;

    uint64_t q = FindFirstPrimeUp(log_q, 2*k);

    cout << "q = " << q << endl;
    
    ZZ a, b, c;
    Polynomial ax(k, q);
    Polynomial bx(k, q);
    Polynomial cx(k, q);


    cout << "ax N = " << ax.GetN() << endl;


    a = RandomBits_ZZ(osize_bit);
    b = RandomBits_ZZ(osize_bit);

    c = a * b;

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "c = " << c << endl;

    TransformZZtoCPoly(ax, a);
    TransformZZtoCPoly(bx, b);

    cout << "ax = " << ax << endl;
    cout << "bx = " << bx << endl;

    ax.SetFormatEval();
    bx.SetFormatEval();

    cx = ax * bx;

    cx.SetFormatCoef();

    cout << "cx = " << cx << endl;

    ZZ crec;
    ReconstructZZfromCPoly(crec, cx);

    cout << "crec = " << crec << endl;
    cout << "crec-c = " << crec-c << endl;
    cout << "c-crec = " << c-crec << endl;

    return 0;
}

void TransformZZtoCPoly(Polynomial &ax, const ZZ a) {
    ZZ a_copy = a;
    for (size_t i = 0; i < ax.GetN(); i++) {
        ax[i] = a_copy % (1<<M);
        a_copy /= (1<<M);
    }
}

void ReconstructZZfromCPoly(ZZ &a, const Polynomial ax) {
    a = 0;
    uint64_t x(1<<M);
    for (int i = ax.GetN() - 1; i >= 0; i--) {
        a = (a * x + ax.at(i));
    }
}