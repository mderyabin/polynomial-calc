#include "examples_utils.h"

#include "polynomial.h"

using namespace std;
using namespace polycalc;

int main(int argc, char const *argv[]) {
        // only independent parameters, specified by user
    size_t logN = 3;
    size_t logq = 5;

    // read parameters from command line (logN first, logq second)
    read_command_line(logN, logq, argc, argv);

    size_t N = 1 << logN;
    size_t M = N << 1;
    uint64_t q = FindFirstPrimeUp(logq, M);

    cout << "logq = " << logq << endl;
    cout << "logN = " << logN << endl;
    cout << "N = " << N << endl;
    cout << "M = " << M << endl;
    cout << "q = " << q << endl;

    Polynomial A(N, q);
    A.GenerateUniform();

    Polynomial B(N, q);

    B = -A;

    cout << "A(x) = " << A << endl;

    cout << "B(x) = " << B << endl;


    cout << "A(x) + B(x) = " << A + B << endl;

    Polynomial C(N, q);
    C.GenerateUniform();

    Polynomial D = A - C;

    cout << "D(x) = " << D << endl;

    cout << "D(X) + C(X) = " << D + C << endl;

    cout << "D(x) = " << D << endl;

    cout << "A(x) = " << A << endl;

    uint64_t c = 5;

    D = c * A;

    cout << c << " * A(x) = " << D << endl;

    uint64_t cinv = ModInvPrime(c, q);

    cout << c << " * A(x) * " << cinv << "  = " << D * cinv << endl;

    return 0;
}