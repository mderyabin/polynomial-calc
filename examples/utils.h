#include <iostream>
#include <string>

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