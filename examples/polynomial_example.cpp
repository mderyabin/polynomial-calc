#include <iostream>

#include "polynomial.h"

using namespace std;

void read_command_line(size_t &logN, size_t &logm, int argc, char const *argv[]);

int main(int argc, char const *argv[])
{

    // only independent parameters, specified by user
    size_t logN = 3;
    size_t logm = 5;

    // read parameters from command line (logN first, logm second)
    read_command_line(logN, logm, argc, argv);

    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t m = FindFirstPrimeUp(logm, M);

    cout << "logm = " << logm << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;
    cout << "m = " << m << endl;

    Polynomial A(N, m);
    A.GenerateUniform();

    Polynomial B(N, m);
    B.GenerateUniform();

    cout << "A(x) = " << A << endl;
    cout << "B(x) = " << B << endl;

    cout << "naive A(x) * B(x) = " << A*B << endl;

    A.SetFormatEval(); // here we do NTT
    B.SetFormatEval(); // here we do NTT

    cout << "NTT(A(x)) = " << A << endl;
    cout << "NTT(B(x)) = " << B << endl;

    Polynomial C = A * B;

    cout << "NTT(C) = NTT(A(x)) * NTT(B(x)) = " << C << endl;

    C.SetFormatCoef(); // here we do inverse NTT

    cout << "C = INTT(NTT(C)) = " << C << endl;
	
    return 0;
}

void read_command_line(size_t &logN, size_t &logm, int argc, char const *argv[]) {
    if (argc == 1) {
        cout << "Parameters logN and logm are set to default." << endl;
        cout << "To set this parameters from command line, use template: " << endl;
        cout << argv[0] << " [logN] [logm]" << endl << endl;
    }

    if (argc > 1) {
        logN = stoi(string(argv[1]));
    }
    if (argc > 2) {
        logm = stoi(string(argv[2]));
    }
}
