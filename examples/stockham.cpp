#include <iostream>

#include "examples_utils.h"
#include "polymath.h"
#include "ntt.h"

using namespace std;
using namespace polycalc;
using namespace NTL;

void printpoly(uint64_t *poly, size_t N, string polyname);
void printvec(uint64_t *vec, size_t N, string vecname);

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

    uint64_t *tf_br = new uint64_t[N];
    uint64_t *pows = new uint64_t[N];


    ComputeTwiddleStockham(tf_br, N, q, false);
    ComputeNWCSequence(pows, N, q, false);

    // ComputeTwiddleFactorsNaive(tf_br, N, q, false);

    printvec(tf_br, N, "tf_br");
    printvec(pows, N, "pows");


    uint64_t *A = new uint64_t[N];
    uint64_t *B = new uint64_t[N];
    GenerateUniformPoly(A, N, q);  
    printpoly(A, N, "A");

    NTTInstance ntt = NTTManager::GetNTTPtr(N, q);

    // ntt->ComputeForward(B, A);
    NaiveNTT(B, A, tf_br, pows, N, q);
    // ntt->ComputeInverse(B);

    printvec(B, N, "B");

    // NWC(A, A, pows, N, q);
    StockhamNTT(A, tf_br, N, q, prec, logq);

    printvec(A, N, "A");

    
    delete [] pows;

    delete [] B;
    delete [] A;
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