#include <iostream>

#include "polymath.h"

using namespace std;

int main(int argc, char const *argv[]) {
    size_t logN = 3;
    size_t logm = 5;
    size_t N = 1 << logN;
    size_t M = N << 1;

    cout << "logm = " << logm << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;


    uint64_t m = FindFirstPrimeUp(logm, M);
    cout << "m = " << m << endl;

    uint64_t g = FindGenerator(m, M);

    cout << "g = " << g << endl;


    uint64_t *tf = new uint64_t[N];
    uint64_t *itf = new uint64_t[N];

    
    ComuteTwiddleFactors(tf, N, m);
    ComuteTwiddleFactors(itf, N, m, true);

    for (size_t i = 0; i < N; i++) {
        cout << "tf[" << i << "] = " << tf[i] << endl;
    }

    for (size_t i = 0; i < N; i++) {
        cout << "itf[" << i << "] = " << itf[i] << endl;
    }

    uint64_t *A = new uint64_t[N];
    uint64_t *nttA = new uint64_t[N];
    uint64_t *B = new uint64_t[N];
    uint64_t *nttB = new uint64_t[N];


    uint64_t *nttC = new uint64_t[N];
    uint64_t *C = new uint64_t[N];
    uint64_t *D = new uint64_t[N];

    
    GenerateUniformPoly(A, N, m);  
    cout << "A = ";
    for (size_t i = 0; i < N-1; i++) {
        cout << A[i] << "*X^" << i << " + ";
    }
    cout << endl;
        
    GenerateUniformPoly(B, N, m);  
    cout << "B = ";
    for (size_t i = 0; i < N; i++) {
        cout << B[i] << "*X^" << i << " + ";
    }
    cout << endl;

    NaiveNTT(nttA, A, tf, N, m);
    cout << "NTT(A) = [";
    for (size_t i = 0; i < N; i++) {
        cout << nttA[i] << ", ";
    }
    cout << "]" << endl;

    NaiveNTT(nttB, B, tf, N, m);
    cout << "NTT(B) = [";
    for (size_t i = 0; i < N; i++) {
        cout << nttB[i] << ", ";
    }
    cout << "]" << endl;

    ModHadamardMul(nttC, nttA, nttB, N, m);

    cout << "NTT(C) = NTT(A * B) = [";
    for (size_t i = 0; i < N; i++) {
        cout << nttC[i] << ", ";
    }
    cout << "]" << endl;

    cout << "C = iNTT(NTT(C)) = ";
    NaiveInvNTT(C, nttC, itf, N, m);
    for (size_t i = 0; i < N; i++) {
        cout << C[i] << "*X^" << i << " + ";
    }
    cout << endl;
    
    NaiveNegacyclicConvolution(D, A, B, N, m);

    cout << "D = naive A*B = ";
    for (size_t i = 0; i < N; i++) {
        cout << D[i] << "*X^" << i << " + ";
    }
    cout << endl;

    delete [] nttC;
    delete [] C;
    delete [] D;
    delete [] nttB;
    delete [] B;
    delete [] nttA;
    delete [] A;
    delete [] itf;
    delete [] tf;
    return 0;
}