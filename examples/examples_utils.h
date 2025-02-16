#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <complex>

#include <NTL/ZZ.h>


long double generaternd(long double a, long double b) {
    std::random_device rd;
    std::uniform_real_distribution<long double> dist(a, b);
    std::default_random_engine en(rd());
    return dist(en);
}

void read_command_line(size_t &logN, size_t &logq, int argc, char const *argv[]) {
    if (argc == 1) {
        std::cout << "Parameters logN and logq are set to default." << std::endl;
        std::cout << "To set this parameters from command line, use template: " << std::endl;
        std::cout << argv[0] << " [logN] [logq]" << std::endl << std::endl;
    }

    if (argc > 1) {
        logN = std::stoi(std::string(argv[1]));
    }
    if (argc > 2) {
        logq = std::stoi(std::string(argv[2]));
    }
}

void TransformZZtoCPoly(uint64_t *ax, const NTL::ZZ a, size_t n, size_t s) {    
    NTL::ZZ a_copy = a;
    for (size_t i = 0; i < n; i++) {
        ax[i] = a_copy % (1<<s);
        a_copy /= (1<<s);
    }
}

void ReconstructZZfromCPoly(NTL::ZZ &a, const uint64_t *ax, size_t n, size_t s) {    
    a = 0;
    uint64_t x(1<<s);
    a = 0;
    for (int i = n - 1; i >= 0; i--) {
        a = (a * x + ax[i]);
    }
}

void PrintArray(const uint64_t *ax, size_t n) {
    std::cout << "[" << ax[0];
    for (size_t i = 1; i < n; i++) {
        std::cout << ", " << ax[i];
    }
    std::cout << "]" << std::endl;
}
