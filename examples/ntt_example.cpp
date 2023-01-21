#include <iostream>
#include <string>

#include "polymath.h"

using namespace std;

void printpoly(uint64_t *poly, size_t N, string polyname);
void printvec(uint64_t *vec, size_t N, string vecname);

void read_command_line(size_t &logN, size_t &logm, int argc, char const *argv[]);

int main(int argc, char const *argv[]) {

    // only independent parameters, specified by user
    size_t logN = 3;
    size_t logm = 5;

    // read parameters from command line (logN first, logm second)
    read_command_line(logN, logm, argc, argv);

    size_t N = 1 << logN;
    size_t M = N << 1;

    cout << "logm = " << logm << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;

    uint64_t m = FindFirstPrimeUp(logm, M);
    cout << "m = " << m << endl;

    uint64_t prec = BarrettPrecompute(m, logm);

    uint64_t invN = ModInvPrime(N, m);
    uint64_t prec_invN = ShoupPrecompute(invN, m);

    uint64_t g = FindGenerator(m, M);

    cout << "g = " << g << endl;

    uint64_t *tf_br = new uint64_t[N];
    uint64_t *itf_br = new uint64_t[N];

    ComputeTwiddleFactors(tf_br, N, m, false);
    ComputeTwiddleFactors(itf_br, N, m, true);

    // for (size_t i = 0; i < N; i++) {
    //     cout << "tf_br[" << i << "] = " << tf_br[i] << endl;
    // }

    // for (size_t i = 0; i < N; i++) {
    //     cout << "itf_br[" << i << "] = " << itf_br[i] << endl;
    // }

    uint64_t *A = new uint64_t[N];
    uint64_t *B = new uint64_t[N];
    uint64_t *C = new uint64_t[N];
    uint64_t *D = new uint64_t[N];

    GenerateUniformPoly(A, N, m);  
    printpoly(A, N, "A");

    GenerateUniformPoly(B, N, m);  
    printpoly(B, N, "B");

    NaiveNegacyclicConvolution(D, A, B, N, m);
    printpoly(D, N, "D = naive A*B");

    CooleyTukeyForwardNTT(A, tf_br, N, m, prec, logm);
    printvec(A, N, "NTT(A)");

    CooleyTukeyForwardNTT(B, tf_br, N, m, prec, logm);
    printvec(B, N, "NTT(B)");

    ModHadamardMul(C, A, B, N, m);
    printvec(C, N, "NTT(C) = NTT(A) * NTT(B)");

    GentlemanSandeInverseNTT(C, itf_br, N, m, invN, prec, prec_invN, logm);
    printpoly(C, N, "C = iNTT(NTT(C))");

    CooleyTukeyForwardNTT(D, tf_br, N, m, prec, logm);
    printvec(D, N, "NTT(D)"); 

    delete [] C;
    delete [] D;
    delete [] B;
    delete [] A;
    delete [] itf_br;
    delete [] tf_br;
    return 0;
}

void printvec(uint64_t *vec, size_t N, string vecname) {
    cout << vecname << " = [";
    cout << vec[0];
    for (size_t i = 1; i < N; i++) {
        cout << ", " << vec[i];
    }
    cout << "]" << endl;
}

void printpoly(uint64_t *poly, size_t N, string polyname) {
    cout << polyname << " = ";
    cout << poly[0];
    for (size_t i = 1; i < N; i++) {
        cout << " + " << poly[i] << "*X^" << i;
    }
    cout << endl;
}

void read_command_line(size_t &logN, size_t &logm, int argc, char const *argv[]) {
    if (argc > 1) {
        logN = stoi(string(argv[1]));
    }
    if (argc > 2) {
        logm = stoi(string(argv[2]));
    }
}