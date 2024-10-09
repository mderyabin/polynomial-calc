// https://dl.acm.org/doi/pdf/10.1145/1277548.1277572

#include <iostream>
#include "encode.h"
#include "ntmath.h"
#include "polynomial.h"

#include <NTL/ZZ.h>
#include <fftw3.h>

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

// const size_t N = 16;

void print_poly(const double *ax, size_t NN);
void print_complex_array(const fftw_complex *a, size_t NN);
void GenerateUniformPoly(double *ax, size_t NN, double lo_b = -10.0, double hi_b = 10.0);
void FFTWPolyMult(double *cx, double *ax, double *bx, size_t NN);
void FFTWPolyMultR(double *cx, double *ax, double *bx, size_t NN);
void NaiveNegacyclicConvolution(double *cx, const double *ax, const double *bx, size_t NN);

void TransformZZtoPoly(double *ax, const ZZ a, size_t NN);
void ReconstructZZfromPoly(ZZ &a, double *ax, size_t NN);

int main() {

    cout << "osize_bit = " << osize_bit << endl;
    cout << "N = " << N << endl;
    cout << "w = " << w << endl;
    cout << "M = " << M << endl;
    cout << "k = " << k << endl;

    ZZ a, b, c;

    a = RandomBits_ZZ(osize_bit);
    b = RandomBits_ZZ(osize_bit);

    c = a * b;

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "c = " << c << endl;

    double *ax = new double[k];
    double *bx = new double[k];
    double *cx = new double[k];

    TransformZZtoPoly(ax, a, k);
    TransformZZtoPoly(bx, b, k);

    cout << endl;
    cout << "ax(X) = "; print_poly(ax, k); cout << endl << endl;
    cout << "bx(X) = "; print_poly(bx, k); cout << endl << endl;

    FFTWPolyMult(cx, ax, bx, k);

    cout << "cx(X) = "; print_poly(cx, k); cout << endl << endl;

    ZZ c_rec; 
    ReconstructZZfromPoly(c_rec, cx, k);

    cout << "c = " << c << endl;
    cout << "c_rec = " << c_rec << endl;

    delete [] cx;
    delete [] bx;
    delete [] ax;

    // double ax[N], bx[N], cx[N], cx_naive[N];
    // GenerateUniformPoly(ax, N);
    // GenerateUniformPoly(bx, N);

    // cout << "a(X) = "; print_poly(ax, N); cout << endl;
    // cout << "b(X) = "; print_poly(bx, N); cout << endl;

    // FFTWPolyMult(cx, ax, bx, N);

    // cout << "c(X) = "; print_poly(cx, N); cout << endl;

    // NaiveNegacyclicConvolution(cx_naive, ax, bx, N);
    // cout << "c_naive(X) = "; print_poly(cx_naive, N); cout << endl;

    return 0;
}

void FFTWPolyMult(double *cx, double *ax, double *bx, size_t NN) {
    fftw_complex *ax_fft, *bx_fft, *cx_fft;
    fftw_plan p, q;

    ax_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NN); 
    bx_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NN);
    cx_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NN);

    p = fftw_plan_dft_r2c_1d(NN, ax, ax_fft, FFTW_ESTIMATE);

    fftw_execute(p);
    fftw_execute_dft_r2c(p, bx, bx_fft);

    cout << "ax_fntt: " << endl;
    print_complex_array(ax_fft, NN);
    cout << endl;

    cout << "bx_fntt: " << endl;
    print_complex_array(bx_fft, NN);
    cout << endl;


    for (size_t i = 1; i < NN/2; i++) {
        cx_fft[i][0] = ax_fft[i][0]*bx_fft[i][0] - ax_fft[i][1]*bx_fft[i][1];
        cx_fft[i][1] = ax_fft[i][0]*bx_fft[i][1] + ax_fft[i][1]*bx_fft[i][0];
    }
    cx_fft[0][0] = ax_fft[0][0]*bx_fft[0][0];
    cx_fft[NN/2][0] = ax_fft[NN/2][0]*bx_fft[NN/2][0];

    cout << "cx_fntt: " << endl;
    print_complex_array(cx_fft, NN);
    cout << endl;

    
    q = fftw_plan_dft_c2r_1d(NN, cx_fft, cx, FFTW_ESTIMATE);
    fftw_execute(q);

    double invN = 1.0 / NN;
    for (size_t i = 0; i < NN; i++) {
        cx[i] *= invN; 
    }

    fftw_destroy_plan(p);
    fftw_destroy_plan(q);

    fftw_free(cx_fft);
    fftw_free(bx_fft);
    fftw_free(ax_fft);

    fftw_cleanup();
}

void GenerateUniformPoly(double *ax, size_t NN, double lo_b, double hi_b) {
    std::random_device rd;
    std::default_random_engine en(rd());
    std::uniform_real_distribution<double> rnd(lo_b, hi_b);

    for (size_t i = 0; i < NN; i++) { 
        ax[i] = rnd(en);
    }
}

void print_poly(const double *ax, size_t NN) {
    cout << ax[0] << " + " << ax[1] << "X"; 
    for (size_t i = 2; i < NN; i++) {
        cout << " + " << ax[i] << "X^" << i; 
    }
}

void NaiveNegacyclicConvolution(double *cx, const double *ax, const double *bx, size_t NN) {
    for (size_t i = 0; i < NN; i++) {
        cx[i] = 0;
        for (size_t j = 0; j <= i; j++) {
            cx[i] += ax[j] * bx[i - j];
        }
        for (size_t j = i+1; j <= N-1; j++) {
            cx[i] += ax[j] * bx[N + i - j];
        }
    }   
}

void print_complex_array(const fftw_complex *a, size_t NN) {
    for (size_t i = 0; i < NN; i++) {
        cout << " " << i << ": " << a[i][0] << " + " << a[i][1] << "*I" << endl;
    }
}

void TransformZZtoPoly(double *ax, const ZZ a, size_t NN) {
    ZZ a_copy = a;
    for (size_t i = 0; i < NN; i++) {
        ax[i] = static_cast<double>(a_copy % (1<<M));
        a_copy /= (1<<M);
    }
}

void ReconstructZZfromPoly(ZZ &a, double *ax, size_t NN) {
    a = 0;
    uint64_t x(1<<M);
    for (int i = NN - 1; i >= 0; i--) {
        a = (a * x + static_cast<uint64_t>(round(ax[i])));
    }
}