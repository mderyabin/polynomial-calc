#include <iostream>
#include <chrono>


#include "polynomial.h"

using namespace std;
using namespace std::chrono;

void read_command_line(size_t &logN, size_t &logm, int argc, char const *argv[]);

void print_time(steady_clock::time_point start, steady_clock::time_point stop, string comment);

int main(int argc, char const *argv[]) {
    steady_clock::time_point start, stop;

    // only independent parameters, specified by user
    size_t logN = 12;
    size_t logq = 58;

    // read parameters from command line (logN first, logm second)
    read_command_line(logN, logq, argc, argv);

    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeDown(logq, M);

    cout << "logq = " << logq << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;
    cout << "q = " << q << endl;

    Polynomial A(N, q);
    Polynomial B(N, q);
    Polynomial C(N, q);

    A.GenerateUniform();
    B.GenerateUniform();

    cout << endl << "NTT transformations performance evaluation" << endl; 

    // segment for Naive NTT  
    uint64_t *tf = new uint64_t[N];
    uint64_t *itf = new uint64_t[N];
    uint64_t *pows = new uint64_t[N];
    uint64_t *ipows = new uint64_t[N];
    uint64_t *ax = new uint64_t[N];
    uint64_t *nttax = new uint64_t[N];
    uint64_t *inttax = new uint64_t[N];

    GenerateUniformPoly(ax, N, q); 

    ComputeTwiddleFactorsNaive(tf, N, q, false);
    ComputeTwiddleFactorsNaive(itf, N, q, true);
    ComputeNWCSequence(pows, N, q, false);
    ComputeNWCSequence(ipows, N, q, true);

    start = high_resolution_clock::now();
    NaiveNTT(nttax, ax, tf, pows, N, q);
    stop = high_resolution_clock::now();
    print_time(start, stop, "Naive NTT");

    start = high_resolution_clock::now();
    NaiveInvNTT(inttax, nttax, itf, ipows, N, q);
    stop = high_resolution_clock::now();
    print_time(start, stop, "Naive Inverse NTT");

    delete [] inttax;
    delete [] nttax;
    delete [] ax;
    delete [] ipows;
    delete [] pows;
    delete [] itf;
    delete [] tf;

    // segment for actual NTT 
    start = high_resolution_clock::now();
    A.SetFormatEval();
    stop = high_resolution_clock::now();
    print_time(start, stop, "Cooley-Tukey Forward NTT");

    start = high_resolution_clock::now();
    A.SetFormatCoef();
    stop = high_resolution_clock::now();
    print_time(start, stop, "Gentleman-Sande Inverse NTT");

    cout << endl << "Polynomial multiplication performance evaluation" << endl; 

    // Segment for Naive multiplication
    start = high_resolution_clock::now();
    C = A * B;
    stop = high_resolution_clock::now();
    print_time(start, stop, "Naive multiplication");

    // Segment for NTT multiplication (includes transforms)
    start = high_resolution_clock::now();
    A.SetFormatEval();
    B.SetFormatEval();
    C = A * B;
    C.SetFormatCoef();
    stop = high_resolution_clock::now();
    print_time(start, stop, "NTT-based multiplication (with transforms)");

    A.SetFormatEval();
    B.SetFormatEval();

    // Segment for NTT multiplication (includes transforms)
    start = high_resolution_clock::now();
    C = A * B;
    stop = high_resolution_clock::now();
    print_time(start, stop, "NTT-based multiplication (without transforms)");

    return 0;
}

void read_command_line(size_t &logN, size_t &logm, int argc, char const *argv[]) {
    if (argc == 1) {
        cout << "Parameters logN and logq are set to default." << endl;
        cout << "To set this parameters from command line, use template: " << endl;
        cout << argv[0] << " [logN] [logq]" << endl << endl;
    }

    if (argc > 1) {
        logN = stoi(string(argv[1]));
    }
    if (argc > 2) {
        logm = stoi(string(argv[2]));
    }
}

void print_time(steady_clock::time_point start, steady_clock::time_point stop, string comment) {
    auto duration = duration_cast<microseconds>(stop - start);
 
    cout << comment << ": " << duration.count() << " microseconds" << endl;
}