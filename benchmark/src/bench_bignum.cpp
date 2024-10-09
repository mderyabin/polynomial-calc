#include <iostream>
#include "encode.h"
#include "ntmath.h"
#include "polynomial.h"

#include "benchmark/benchmark.h"


#include <NTL/ZZ.h>
#include <fftw3.h>


using namespace std;
using namespace polycalc;
using namespace NTL;

// const size_t osize_bit = 1024;

// const size_t N = 2*osize_bit;       // size of numbers (product of a and b size) 
// const size_t w = 64;                // word size

// const size_t M = 16;                 // coefficients size
// const size_t k = N/M;               // number of coefficients

// const size_t log_q = 58;          

void TransformZZtoCPoly(Polynomial &ax, const ZZ a, size_t M) {
    ZZ a_copy = a;
    for (size_t i = 0; i < ax.GetN(); i++) {
        ax[i] = a_copy % (1<<M);
        a_copy /= (1<<M);
    }
}

void ReconstructZZfromCPoly(ZZ &a, const Polynomial ax, size_t M) {
    a = 0;
    uint64_t x(1<<M);
    for (int i = ax.GetN() - 1; i >= 0; i--) {
        a = (a * x + ax.at(i));
    }
}

void TransformZZtoPoly(double *ax, const ZZ a, size_t NN, size_t M) {
    ZZ a_copy = a;
    for (size_t i = 0; i < NN; i++) {
        ax[i] = static_cast<double>(a_copy % (1<<M));
        a_copy /= (1<<M);
    }
}

void ReconstructZZfromPoly(ZZ &a, double *ax, size_t NN, size_t M) {
    a = 0;
    uint64_t x(1<<M);
    for (int i = NN - 1; i >= 0; i--) {
        a = (a * x + static_cast<uint64_t>(round(ax[i])));
    }
}

void BC_ZZBigNumMul(benchmark::State& state, size_t N) {
    ZZ a, b, c;
    a = RandomBits_ZZ(N/2);
    b = RandomBits_ZZ(N/2);

    while (state.KeepRunning()) {
        c = a * b;
    }
}

void BC_PolyBigNumMul(benchmark::State& state, size_t N, size_t M, size_t log_q) {
    size_t k = N/M;
    uint64_t q = FindFirstPrimeUp(log_q, 2*k);
    Polynomial ax(k, q);
    Polynomial bx(k, q);
    Polynomial cx(k, q);

    ZZ a, b, c;
    a = RandomBits_ZZ(N/2);
    b = RandomBits_ZZ(N/2);

    TransformZZtoCPoly(ax, a, M);
    TransformZZtoCPoly(bx, b, M);

    while (state.KeepRunning()) {
        ax.SetFormatEval();
        bx.SetFormatEval();
        cx = ax * bx;
        cx.SetFormatCoef();
    }
}

void BC_PolyBigNumMulFFTW(benchmark::State& state, size_t N, size_t M) {
    size_t k = N/M;

    double *ax = new double[k];
    double *bx = new double[k];
    double *cx = new double[k];

    fftw_complex ax_fft[k], bx_fft[k], cx_fft[k];
    fftw_plan p, q;

    p = fftw_plan_dft_r2c_1d(k, ax, ax_fft, FFTW_ESTIMATE);
    q = fftw_plan_dft_c2r_1d(k, cx_fft, cx, FFTW_ESTIMATE);
    double invN = 1.0 / k;



    // ax_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * k); 
    // bx_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * k);
    // cx_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * k);


    ZZ a, b, c;
    a = RandomBits_ZZ(N/2);
    b = RandomBits_ZZ(N/2);

    TransformZZtoPoly(ax, a, k, M);
    TransformZZtoPoly(bx, b, k, M);

    while (state.KeepRunning()) {
        fftw_execute(p);
        fftw_execute_dft_r2c(p, bx, bx_fft);
        for (size_t i = 1; i < k/2; i++) {
            cx_fft[i][0] = ax_fft[i][0]*bx_fft[i][0] - ax_fft[i][1]*bx_fft[i][1];
            cx_fft[i][1] = ax_fft[i][0]*bx_fft[i][1] + ax_fft[i][1]*bx_fft[i][0];
        }
        cx_fft[0][0] = ax_fft[0][0]*bx_fft[0][0];
        cx_fft[k/2][0] = ax_fft[k/2][0]*bx_fft[k/2][0];
        fftw_execute(q);
        for (size_t i = 0; i < k; i++) {
            cx[i] *= invN; 
        }
    }

    fftw_destroy_plan(p);
    fftw_destroy_plan(q);

    // fftw_free(cx_fft);
    // fftw_free(bx_fft);
    // fftw_free(ax_fft);

    fftw_cleanup();

    delete [] cx;
    delete [] bx;
    delete [] ax;
}


BENCHMARK_CAPTURE(BC_ZZBigNumMul, BC_ZZBigNumMul_256, 256);
BENCHMARK_CAPTURE(BC_ZZBigNumMul, BC_ZZBigNumMul_512, 512);
BENCHMARK_CAPTURE(BC_ZZBigNumMul, BC_ZZBigNumMul_1024, 1024);

BENCHMARK_CAPTURE(BC_PolyBigNumMul, BC_PolyBigNumMul_256_16_58, 256, 16, 58);
BENCHMARK_CAPTURE(BC_PolyBigNumMul, BC_PolyBigNumMul_512_16_58, 512, 16, 58);
BENCHMARK_CAPTURE(BC_PolyBigNumMul, BC_PolyBigNumMul_1024_16_58, 1024, 16, 58);

BENCHMARK_CAPTURE(BC_PolyBigNumMul, BC_PolyBigNumMul_256_8_58, 256, 8, 58);
BENCHMARK_CAPTURE(BC_PolyBigNumMul, BC_PolyBigNumMul_512_8_58, 512, 8, 58);
BENCHMARK_CAPTURE(BC_PolyBigNumMul, BC_PolyBigNumMul_1024_8_58, 1024, 8, 58);

BENCHMARK_CAPTURE(BC_PolyBigNumMulFFTW, BC_PolyBigNumMulFFTW_256_32, 256, 32);
BENCHMARK_CAPTURE(BC_PolyBigNumMulFFTW, BC_PolyBigNumMulFFTW_512_32, 512, 32);
BENCHMARK_CAPTURE(BC_PolyBigNumMulFFTW, BC_PolyBigNumMulFFTW_1024_32, 1024, 32);


BENCHMARK_CAPTURE(BC_PolyBigNumMulFFTW, BC_PolyBigNumMulFFTW_256_16, 256, 16);
BENCHMARK_CAPTURE(BC_PolyBigNumMulFFTW, BC_PolyBigNumMulFFTW_512_16, 512, 16);
BENCHMARK_CAPTURE(BC_PolyBigNumMulFFTW, BC_PolyBigNumMulFFTW_1024_16, 1024, 16);

BENCHMARK_CAPTURE(BC_PolyBigNumMulFFTW, BC_PolyBigNumMulFFTW_256_8, 256, 8);
BENCHMARK_CAPTURE(BC_PolyBigNumMulFFTW, BC_PolyBigNumMulFFTW_512_8, 512, 8);
BENCHMARK_CAPTURE(BC_PolyBigNumMulFFTW, BC_PolyBigNumMulFFTW_1024_8, 1024, 8);

BENCHMARK_MAIN();
