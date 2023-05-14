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
    cout << "A(x) = " << A << endl;

    A.Serialize("poly.txt", JSON);

    Polynomial B;
    B.Deserialize("poly.txt", JSON);

    cout << "B(x) = " << B << endl;

    cout << "N = " << B.GetN() << endl;
    cout << "q = " << B.GetModulus() << endl;
}