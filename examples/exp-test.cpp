#include <iostream>
#include <string>
#include <cmath>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

ZZ FastModExp(ZZ a, ZZ e, ZZ m);

int main(int argc, char const *argv[]) {
    size_t m_bsize = 128;

    if (argc > 1) {
        int r = stoi(argv[1]); 
        m_bsize = size_t(1) << static_cast<size_t>(ceil(log2(double(r))));
    }

    ZZ m;
    ZZ a, e, c, inv;

    GenPrime(m, m_bsize);

    cout << "m_bsize = " << m_bsize << endl;
    cout << "m = " << m << endl;

    a = RandomBits_ZZ(m_bsize - 1);
    e = RandomBits_ZZ(m_bsize - 1);

    cout << "a = " << a << endl;
    cout << "e = " << e << endl;

    c = PowerMod(a, e, m);

    cout << "a^e mod m = " << c << endl;

    inv = PowerMod(a, m-2, m);


    cout << "a^(m-2) mod m = " << inv << endl;
    cout << "a*c mod m = " << MulMod(a, inv, m) << endl;

    ZZ c1, c2;

    c1 = FastModExp(a, e, m);

    cout << " -- testing -- " << endl;
    cout << "a^e mod m = " << c1 << endl;

    c2 = FastModExp(a, m-2, m);

    cout << "a^(m-2) mod m = " << c2 << endl;
    cout << "a*c mod m = " << MulMod(a, c2, m) << endl;

    cout << " -   diff c: " << c1 - c << endl;
    cout << " - diff inv: " << c2 - inv << endl;

    return 0;
}

ZZ FastModExp(ZZ a, ZZ e, ZZ m) { 
    ZZ res, base, ee;
    bool flag = false;
    base = a;
    ee = e;

    while (ee > 0) {
        if (ee % 2 == 1) {
            if (flag) { 
                res = MulMod(res, base, m);
            } else {
                res = base;
                flag = true;
            }
        }
        ee >>= 1;
        if (ee != 0) {
            base = SqrMod(base, m);
        }
    }

    return res;
}