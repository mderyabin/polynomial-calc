#include <iostream>
#include <string>

#include "vectbigint.h"
#include "examples_utils.h"

using namespace std;
using namespace NTL;
using namespace polycalc;


ZZ FastModExp(ZZ a, ZZ e, ZZ m);

// a >= b ??

int main(int argc, char const *argv[]) {
    size_t bsize = 128;
    size_t shift = 1;
    size_t s = 16;

    if (argc > 1) {
        int r = stoi(argv[1]); 
        bsize = size_t(1) << static_cast<size_t>(ceil(log2(double(r))));
    }

    if (argc > 2) {
        shift = stoi(argv[2]); 
    }

    if (argc > 3) {
        s = stoi(argv[3]); 
    }

    size_t n = bsize / s;
    size_t m_bitsize = bsize-s;

    ZZ a, b, m, c; 
    uint64_t *ax = new uint64_t[n];
    uint64_t *bx = new uint64_t[n];
    uint64_t *cx = new uint64_t[n];
    uint64_t *mx = new uint64_t[n];

    
    // globalMemPool.Init(n);
    
    m = RandomBits_ZZ(m_bitsize);
    a = RandomBits_ZZ(m_bitsize) % m;
    b = RandomBits_ZZ(m_bitsize) % m;
    while (b >= a) {
        b = RandomBits_ZZ(m_bitsize) % m;
    }

    
    cout << "bsize = " << bsize << endl;
    cout << "s = " << s << endl;
    cout << "n = " << n << endl;
    cout << "shift = " << shift << endl;
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "m = " << m << endl;

    
    TransformZZtoCPoly(ax, a, n, s);
    TransformZZtoCPoly(bx, b, n, s);
    TransformZZtoCPoly(mx, m, n, s);

    cout << "ax = "; PrintArray(ax, n);
    cout << "bx = "; PrintArray(bx, n);
    cout << "mx = "; PrintArray(mx, n);

    LeftShift(cx, ax, shift, n, s);
    PrintArray(cx, n);
    ReconstructZZfromCPoly(c, cx, n, s);

    cout << " -- shift test -- " << endl;
    cout << "a << " << shift << " = " << (a<<shift) << endl;
    cout << "a << " << shift << " = " << (a<<shift)%(ZZ(1)<<bsize) << endl;
    cout << "cx = "; PrintArray(cx, n);
    cout << "c = " << c << endl;


    copy(ax, ax+n, cx);
    LeftShift(cx, ax, shift, n, s);
    PrintArray(cx, n);
    ReconstructZZfromCPoly(c, cx, n, s);

    cout << " -- shift test inlace -- " << endl;
    cout << "a << " << shift << " = " << (a<<shift) << endl;
    cout << "a << " << shift << " = " << (a<<shift)%(ZZ(1)<<bsize) << endl;
    cout << "cx = "; PrintArray(cx, n);
    cout << "c = " << c << endl;
   

    cout << " -- add lazy test -- " << endl;
    AddLazy(cx, ax, bx, n);
    cout << "cx = ax + bx = "; PrintArray(cx, n);

    ReconstructZZfromCPoly(c, cx, n, s);
    cout << "c = a + b = " << c << endl;
    cout << "    a + b = " << a+b << endl;

    CarryPropagation(cx, n, s);
    cout << " -- carry propagation test -- " << endl;
    cout << "cx = ax + bx = "; PrintArray(cx, n);
    ReconstructZZfromCPoly(c, cx, n, s);
    cout << "c = a + b = " << c << endl;


    cout << " -- sub test -- " << endl;
    SubLazy(cx, ax, bx, n, s);
    cout << "cx = ax - bx = "; PrintArray(cx, n);
    ReconstructZZfromCPoly(c, cx, n, s);
    cout << "c = a - b = " << c << endl;
    cout << "    a - b = " << a-b << endl;

    cout << " -- sub test inplace -- " << endl;
    copy(ax, ax+n, cx);
    SubLazy(cx, bx, n, s);
    cout << "cx = ax - bx = "; PrintArray(cx, n);
    ReconstructZZfromCPoly(c, cx, n, s);
    cout << "c = a - b = " << c << endl;
    cout << "    a - b = " << a-b << endl;


    cout << " -- is bigger or equal test -- " << endl;
    bool ageb = IsBiggerOrEqual(ax, bx, n);
    cout << "test 1 IsBiggerOrEqual(ax, bx, n) = " << (ageb ? "ax >= bx" : "ax < bx") << endl;;
    bool bgea = IsBiggerOrEqual(bx, ax, n);
    cout << "test 2 IsBiggerOrEqual(bx, ax, n) = " << (bgea ? "bx >= ax" : "bx < ax") << endl;;

    cout << " -- add mod test -- " << endl;

    AddMod(cx, ax, bx, mx, n, s);
    cout << "cx = ax + bx mod mx = "; PrintArray(cx, n);
    cout << "              a + b = " << (a+b) << endl;
    cout << "          a + b - m = " << (a+b-m) << endl;
    cout << "        a + b mod m = " << (a+b)%m << endl;
    ReconstructZZfromCPoly(c, cx, n, s);
    cout << "    c = a + b mod m = " << c << endl;

    cout << " -- add mod test inplace -- " << endl;

    copy(ax, ax+n, cx);
    AddMod(cx, bx, mx, n, s);
    cout << "cx = ax + bx mod mx = "; PrintArray(cx, n);
    cout << "              a + b = " << (a+b) << endl;
    cout << "          a + b - m = " << (a+b-m) << endl;
    cout << "        a + b mod m = " << (a+b)%m << endl;
    ReconstructZZfromCPoly(c, cx, n, s);
    cout << "    c = a + b mod m = " << c << endl;

    ZZ tmpp = (a+b)%m;
    TransformZZtoCPoly(cx, tmpp, n, s);
    cout << "(a+b)%m = "; PrintArray(cx, n);

    cout << " -- mul mod test -- " << endl;
    MulMod(cx, ax, bx, mx, n, s);
    cout << "cx = ax + bx mod mx = "; PrintArray(cx, n);
    cout << "        a * b mod m = " << MulMod(a, b, m) << endl;
    ReconstructZZfromCPoly(c, cx, n, s);
    cout << "    c = a * b mod m = " << c << endl;
    cout << "               diff = " << c-MulMod(a, b, m) << endl;

    tmpp = MulMod(a, b, m);
    TransformZZtoCPoly(cx, tmpp, n, s);
    cout << "(a*b)%m = "; PrintArray(cx, n);

    cout << " -- exp mod test -- " << endl;
    FastModExp(cx, ax, bx, mx, n, s);
    cout << "cx = ax ^ bx mod mx = "; PrintArray(cx, n);
    cout << "        a ^ b mod m = " << FastModExp(a, b, m) << endl;
    ReconstructZZfromCPoly(c, cx, n, s);
    cout << "    c = a ^ b mod m = " << c << endl;
    cout << "               diff = " << c-FastModExp(a, b, m) << endl;


    delete [] mx;
    delete [] cx;
    delete [] bx;
    delete [] ax;
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
