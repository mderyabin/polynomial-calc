#include <iostream>
#include <string>

#include "polymath.h"

using namespace std;

void printpoly(uint64_t *poly, size_t N, string polyname);
void printvec(uint64_t *vec, size_t N, string vecname);

void read_command_line(size_t &logN, size_t &logq, int argc, char const *argv[]);

int main(int argc, char const *argv[]) {

    // only independent parameters, specified by user
    size_t logN = 3;
    size_t logq = 5;

    // read parameters from command line (logN first, logq second)
    read_command_line(logN, logq, argc, argv);

    size_t N = 1 << logN;
    size_t M = N << 1;

    cout << "logq = " << logq << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;

    uint64_t q = FindFirstPrimeUp(logq, M);
    cout << "q = " << q << endl;

    uint64_t prec = BarrettPrecompute(q, logq);

    uint64_t invN = ModInvPrime(N, q);
    uint64_t prec_invN = ShoupPrecompute(invN, q);

    uint64_t g = FindGenerator(q, M);

    cout << "g = " << g << endl;

    uint64_t *tf_br = new uint64_t[N];
    uint64_t *itf_br = new uint64_t[N];

    ComputeTwiddleFactors(tf_br, N, q, false);
    ComputeTwiddleFactors(itf_br, N, q, true);

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

    GenerateUniformPoly(A, N, q);  
    printpoly(A, N, "A");

    GenerateUniformPoly(B, N, q);  
    printpoly(B, N, "B");

    NaiveNegacyclicConvolution(D, A, B, N, q);
    printpoly(D, N, "D = naive A*B");

    CooleyTukeyForwardNTT(A, tf_br, N, q, prec, logq);
    printvec(A, N, "NTT(A)");

    CooleyTukeyForwardNTT(B, tf_br, N, q, prec, logq);
    printvec(B, N, "NTT(B)");

    ModHadamardMul(C, A, B, N, q);
    printvec(C, N, "NTT(C) = NTT(A) * NTT(B)");

    GentlemanSandeInverseNTT(C, itf_br, N, q, invN, prec, prec_invN, logq);
    printpoly(C, N, "C = iNTT(NTT(C))");

    CooleyTukeyForwardNTT(D, tf_br, N, q, prec, logq);
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

void read_command_line(size_t &logN, size_t &logq, int argc, char const *argv[]) {
    if (argc == 1) {
        cout << "Parameters logN and logq are set to default." << endl;
        cout << "To set this parameters from command line, use template: " << endl;
        cout << argv[0] << " [logN] [logq]" << endl << endl;
    }

    if (argc > 1) {
        logN = stoi(string(argv[1]));
    }
    if (argc > 2) {
        logq = stoi(string(argv[2]));
    }
}