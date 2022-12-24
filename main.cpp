#include <iostream>

#include "src/include/math.h"

using namespace std;

int main(int argc, char const *argv[])
{
    uint64_t m = 11;
    uint64_t a = 5;
    uint64_t b = 9;

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "m = " << m << endl;

    cout << "a + b mod m = " << ModAdd(a, b, m) << endl;
    
    return 0;
}
