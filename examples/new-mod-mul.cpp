#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

class MemoryPool {
public:
    MemoryPool() {}
    ~MemoryPool();

    void Init(size_t n);

    uint64_t *base;
    uint64_t *shifted;
    uint64_t *tmp1;
    uint64_t *tmp2;
    uint64_t *tmp3;

} globalMemPool;

void TransformZZtoCPoly(uint64_t *ax, const ZZ a, size_t n, size_t s);
void ReconstructZZfromCPoly(ZZ &a, const uint64_t *ax, size_t n, size_t s);
void PrintArray(const uint64_t *ax, size_t n);
void LeftShift(uint64_t *cx, const uint64_t *ax, size_t shift, size_t n, size_t s);
// void LeftShift(uint64_t *cx, size_t shift, size_t n, size_t s);

bool IsBiggerOrEqual(const uint64_t *ax, const uint64_t *bx, size_t n);
uint64_t CarryPropagation(uint64_t *cx, const uint64_t *ax, size_t n, size_t s);
uint64_t CarryPropagation(uint64_t *cx, size_t n, size_t s);
void AddLazy(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t n);
void AddLazy(uint64_t *cx, const uint64_t *ax, size_t n);
uint64_t AddFull(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t n, size_t s);
uint64_t AddFull(uint64_t *cx, const uint64_t *ax, size_t n, size_t s);
uint64_t SubLazy(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t n, size_t s);
uint64_t SubLazy(uint64_t *cx, const uint64_t *bx, size_t n, size_t s);

void AddMod(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s);
void AddMod(uint64_t *cx, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s);

void MulMod(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s);
void MulMod(uint64_t *cx, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s);

void SqrMod(uint64_t *cx, const uint64_t *ax, const uint64_t *mx, size_t n, size_t s);
void SqrMod(uint64_t *cx, const uint64_t *mx, size_t n, size_t s);

void FastModExp(uint64_t *cx, const uint64_t *ax, const uint64_t *ex, const uint64_t *mx, size_t n, size_t s);
ZZ FastModExp(ZZ a, ZZ e, ZZ m);

// a >= b ??
bool IsBiggerOrEqual(const uint64_t *ax, const uint64_t *bx, size_t n);

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

    globalMemPool.Init(n);
    
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

void LeftShift(uint64_t *cx, const uint64_t *ax, size_t shift, size_t n, size_t s) {
    size_t bs = shift % s;
    size_t bc = shift / s;

    // cout << "block_count = " << bc << endl;
    // cout << "bit_shift = " << bs << endl;
    
    size_t mask = (1 << s) - 1;
    uint64_t clo, chi;

    cx[bc] = (ax[0] << bs) & mask;
    for (size_t i = 1; i < n-bc; i++) {
        clo = (ax[i] << bs) & mask;
        chi = ax[i-1] >> (s - bs);
        cx[bc + i] = clo + chi;
    }
}

void LeftShift(uint64_t *cx, size_t shift, size_t n, size_t s) {

    // uint64_t *ax = new uint64_t[n];
    uint64_t *ax = globalMemPool.tmp2;

    size_t bs = shift % s;
    size_t bc = shift / s;

    copy(cx, cx+n, ax);

    // cout << "block_count = " << bc << endl;
    // cout << "bit_shift = " << bs << endl;
    
    size_t mask = (1 << s) - 1;
    uint64_t clo, chi;

    cx[bc] = (ax[0] << bs) & mask;
    for (size_t i = 1; i < n-bc; i++) {
        clo = (ax[i] << bs) & mask;
        chi = ax[i-1] >> (s - bs);
        cx[bc + i] = clo + chi;
    }

    // delete [] ax;
}



void TransformZZtoCPoly(uint64_t *ax, const ZZ a, size_t n, size_t s) {    
    ZZ a_copy = a;
    for (size_t i = 0; i < n; i++) {
        ax[i] = a_copy % (1<<s);
        a_copy /= (1<<s);
    }
}

void ReconstructZZfromCPoly(ZZ &a, const uint64_t *ax, size_t n, size_t s) {    
    a = 0;
    uint64_t x(1<<s);
    a = 0;
    for (int i = n - 1; i >= 0; i--) {
        a = (a * x + ax[i]);
    }
}

void PrintArray(const uint64_t *ax, size_t n) {
    cout << "[" << ax[0];
    for (size_t i = 1; i < n; i++) {
        cout << ", " << ax[i];
    }
    cout << "]" << endl;
}

bool IsBiggerOrEqual(const uint64_t *ax, const uint64_t *bx, size_t n) {
    size_t i = n-1;
    while (ax[i] == bx[i] and (i < n)) i--;
    if (ax[i] < bx[i]) return false;
    else return true;
}

bool IsBigger(const uint64_t *ax, const uint64_t *bx, size_t n) {
    size_t i = n-1;
    while (ax[i] == bx[i] and (i < n)) i--;
    if (ax[i] > bx[i]) return true;
    else return false;
}


uint64_t CarryPropagation(uint64_t *cx, const uint64_t *ax, size_t n, size_t s) {
    uint64_t carry = 0, cxi;
    size_t mask = (1 << s) - 1;
    for (size_t i = 0; i < n; i++) {
        cxi = ax[i] + carry;
        cx[i] = cxi & mask;
        carry = cxi >> s;
    }
    return carry;
}

uint64_t CarryPropagation(uint64_t *cx, size_t n, size_t s) {
    uint64_t carry = 0, cxi;
    size_t mask = (1 << s) - 1;
    for (size_t i = 0; i < n; i++) {
        cxi = cx[i] + carry;
        cx[i] = cxi & mask;
        carry = cxi >> s;
    }
    return carry;
}

void AddLazy(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t n) {
    for (size_t i = 0; i < n; i++) {
        cx[i] = ax[i] + bx[i];
    }
}

void AddLazy(uint64_t *cx, const uint64_t *ax, size_t n) {
    for (size_t i = 0; i < n; i++) {
        cx[i] = cx[i] + ax[i];
    }
}

uint64_t AddFull(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t n, size_t s) {
    AddLazy(cx, ax, bx, n);
    return CarryPropagation(cx, n, s);
}

uint64_t AddFull(uint64_t *cx, const uint64_t *ax, size_t n, size_t s) {
    AddLazy(cx, ax, n);
    return CarryPropagation(cx, n, s);
}

// a > b!!
uint64_t SubLazy(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t n, size_t s) {
    uint64_t neg_carry = 0;
    uint64_t axi, bxi;
    uint64_t base = 1 << s;
    for (size_t i = 0; i < n; i++) {
        axi = ax[i];
        bxi = bx[i]+neg_carry;
        if (axi < bxi) {
            axi += base;
            neg_carry = 1;
        } else {
            neg_carry = 0;
        }
        cx[i] = axi - bxi;
    }
    return neg_carry;
}

uint64_t SubLazy(uint64_t *cx, const uint64_t *bx, size_t n, size_t s) {
    uint64_t neg_carry = 0;
    uint64_t cxi, bxi;
    uint64_t base = 1 << s;
    for (size_t i = 0; i < n; i++) {
        cxi = cx[i];
        bxi = bx[i] + neg_carry;
        if (cxi < bxi) {
            cxi += base;
            neg_carry = 1;
        } else {
            neg_carry = 0;
        }
        cx[i] = cxi - bxi;
    }
    return neg_carry;
}

void AddMod(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s) {
    AddFull(cx, ax, bx, n, s);
    if (IsBiggerOrEqual(cx, mx, n)) {
        SubLazy(cx, mx, n, s);
    }
}

void AddMod(uint64_t *cx, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s) {
    AddFull(cx, bx, n, s);
    if (IsBiggerOrEqual(cx, mx, n)) {
        SubLazy(cx, mx, n, s);
    }
}

void MulMod(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s) {
    // uint64_t *shifted = new uint64_t[n];
    uint64_t *shifted = globalMemPool.shifted;

    copy(ax, ax+n, shifted);

    if (bx[0] % 2 == 1) {
        copy(ax, ax+n, cx);
    } else {
        for (size_t i = 0; i < n; i++) {
            cx[i] = 0;
        }
    }
   
    for (size_t i = 1; i < n*s; i++) {
        LeftShift(shifted, 1, n, s);

        while (IsBiggerOrEqual(shifted, mx, n)) {
            SubLazy(shifted, mx, n, s);
        }
        
        size_t j = i/s;
        size_t k = i%s;

        if ((bx[j] >> k)%2 == 1) AddMod(cx, shifted, mx, n, s);

    }

    // delete [] shifted;
}

void MulMod(uint64_t *cx, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s) {
    // uint64_t *shifted = new uint64_t[n];

    uint64_t *shifted = globalMemPool.shifted;

    copy(cx, cx+n, shifted);

    if (bx[0] % 2 == 0) {
        for (size_t i = 0; i < n; i++) {
            cx[i] = 0;
        }
    }
   
    for (size_t i = 1; i < n*s; i++) {
        LeftShift(shifted, 1, n, s);

        while (IsBiggerOrEqual(shifted, mx, n)) {
            SubLazy(shifted, mx, n, s);
        }
        
        size_t j = i/s;
        size_t k = i%s;

        if ((bx[j] >> k)%2 == 1) AddMod(cx, shifted, mx, n, s);

    }

    // delete [] shifted;
}

void SqrMod(uint64_t *cx, const uint64_t *ax, const uint64_t *mx, size_t n, size_t s) {
    // copy(ax, ax+n, cx);
    MulMod(cx, ax, ax, mx, n, s);
}

void SqrMod(uint64_t *cx, const uint64_t *mx, size_t n, size_t s) {
    // uint64_t *shifted = new uint64_t[n];
    // uint64_t *bx = new uint64_t[n];

    uint64_t *shifted = globalMemPool.shifted;
    uint64_t *bx = globalMemPool.tmp1;

    copy(cx, cx+n, shifted);
    copy(cx, cx+n, bx);

    if (bx[0] % 2 == 0) {
        for (size_t i = 0; i < n; i++) {
            cx[i] = 0;
        }
    }
   
    for (size_t i = 1; i < n*s; i++) {
        LeftShift(shifted, 1, n, s);

        while (IsBiggerOrEqual(shifted, mx, n)) {
            SubLazy(shifted, mx, n, s);
        }
        
        size_t j = i/s;
        size_t k = i%s;

        if ((bx[j] >> k)%2 == 1) AddMod(cx, shifted, mx, n, s);

    }

    // delete [] bx;
    // delete [] shifted;
}

void FastModExp(uint64_t *cx, const uint64_t *ax, const uint64_t *ex, const uint64_t *mx, size_t n, size_t s) { 

    // uint64_t *base = new uint64_t[n];
    uint64_t *base = globalMemPool.base;

    copy(ax, ax+n, base);

    bool flag = false;

    size_t emsb = 0;
    for (size_t i = 0; i < n; i++) {
        if (ex[i] != 0) emsb = i;
    }
    emsb++;

    for (size_t i = 0; i < emsb*s; i++) {
        size_t j = i/s;
        size_t k = i%s;


        if ((ex[j] >> k)%2 == 1) {
            if (flag) { 
                MulMod(cx, base, mx, n, s);
            } else {
                copy(base, base+n, cx);
                flag = true;
            }
        }

        if (! (j >= emsb and (ex[j] >> k) == 0) ) {
            SqrMod(base, mx, n, s);
        }
    }

    // delete [] base;
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

void MemoryPool::Init(size_t n) {
    base = new uint64_t[n];
    shifted = new uint64_t[n];
    tmp1 = new uint64_t[n];
    tmp2 = new uint64_t[n];
    tmp3 = new uint64_t[n];
}

MemoryPool::~MemoryPool() {
    delete [] tmp3;
    delete [] tmp2;
    delete [] tmp1;
    delete [] shifted;
    delete [] base;
}