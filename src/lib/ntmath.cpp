#include "ntmath.h"

#include <iostream>
#include <random>
using namespace std;

typedef unsigned __int128 uint128_t;

uint64_t MSB(uint64_t t) {
    uint64_t msb = 0;
    while (t > 0) {
        msb++;
        t >>= 1;
    }
    return msb;
}

uint64_t ModAdd(const uint64_t a, const uint64_t b, const uint64_t m) {
    uint64_t c = a + b;
    return (c >= m) ? c - m : c;
}

void ModAddEq(uint64_t &a, const uint64_t b, const uint64_t m) {
    a += b;
    a = (a >= m) ? a - m : a;
}

uint64_t ModSub(const uint64_t a, const uint64_t b, const uint64_t m) {
    return (a > b) ? a - b : m + a - b;
}

void ModSubEq(uint64_t &a, const uint64_t b, const uint64_t m) {
    a = (a > b) ? a - b : m + a - b;
}

//*
// https://core.ac.uk/download/pdf/287482281.pdf
uint64_t BarrettPrecompute(uint64_t m, uint64_t logm) {
    uint64_t bar = 2 * logm + 3;
    return (static_cast<uint128_t>(1) << bar) / m;
}

uint64_t ModMultBarrett(uint64_t a, uint64_t b, uint64_t m, uint64_t prec, uint64_t logm) {
    uint64_t res;
    //int64_t bar = 2 * logm;

    uint128_t mul = static_cast<uint128_t>(a) * b;

    uint128_t tmp1 = mul;
    uint128_t tmp2 = tmp1 >> (logm-2);

    tmp1 = tmp2 * prec;
    tmp2 = tmp1 >> (logm + 5);
    tmp1 = tmp2 * m;

    res = static_cast<uint64_t>(mul - tmp1);

    while (res >= m) res -= m;

    return res;
}


void ModMultBarrettEq(uint64_t &a, uint64_t b, uint64_t m, uint64_t prec, uint64_t logm) {

    uint128_t mul = static_cast<uint128_t>(a) * b;

    uint128_t tmp1 = mul;
    uint128_t tmp2 = tmp1 >> (logm-2);

    tmp1 = tmp2 * prec;
    tmp2 = tmp1 >> (logm + 5);
    tmp1 = tmp2 * m;

    a = static_cast<uint64_t>(mul - tmp1);

    while (a >= m) a -= m;
}

//*/

/* то же самое, только 64-бит версия из HEAAN 
// https://github.com/snucrypto/HEAAN/blob/master/HEAAN/src/RingMultiplier.cpp

// logm > 32 !!
uint64_t BarrettPrecompute(uint64_t m, uint64_t logm) {
    uint64_t bar = 2 * logm;
    return (static_cast<uint128_t>(1) << bar) / m;
}

// logm > 32 !!
uint64_t ModMultBarrett(uint64_t a, uint64_t b, uint64_t m, uint64_t prec, uint64_t logm) {
	uint64_t res;
    uint64_t bar = 2 * logm;
    
    uint128_t mul = static_cast<uint128_t>(a) * b;
	uint64_t abot = static_cast<uint64_t>(mul);
	uint64_t atop = static_cast<uint64_t>(mul >> 64);
	uint128_t tmp = static_cast<uint128_t>(abot) * prec;
	tmp >>= 64;
	tmp += static_cast<uint128_t>(atop) * prec;
	tmp >>= bar - 64;
	tmp *= m;
	tmp = mul - tmp;
	res = static_cast<uint64_t>(tmp);  // mod 2^64
	if (res >= m) res -= m;

    return res;
}

//*/


// метод Шупа применяется для умножения на константу
// https://pdfs.semanticscholar.org/e000/fa109f1b2a6a3e52e04462bac4b7d58140c9.pdf

uint64_t ShoupPrecompute(uint64_t c, uint64_t m) {
    uint128_t w = static_cast<uint128_t>(c);
    w <<= 64;
    w /= m;
    return static_cast<uint64_t>(w);
}

uint64_t ModMulShoup(uint64_t a, uint64_t c, uint64_t m, uint64_t prec) {
    uint64_t res;
   
    uint128_t aa  = static_cast<uint128_t>(a);
    uint128_t mul = aa * c;
    uint128_t tmp = ((aa * prec) >> 64) * m;

    res = static_cast<uint64_t>(mul - tmp);  // mod 2^64

    if (res >= m) res -= m;

    return res;
}

void ModMulShoupEq(uint64_t &a, uint64_t c, uint64_t m, uint64_t prec) {
   
    uint128_t aa  = static_cast<uint128_t>(a);
    uint128_t mul = aa * c;
    uint128_t tmp = ((aa * prec) >> 64) * m;

    a = static_cast<uint64_t>(mul - tmp);  // mod 2^64

    if (a >= m) a -= m;
}


uint64_t ModExp(uint64_t a, uint64_t e, uint64_t m, uint64_t prec, uint64_t logm) {
    if (e == 0) return 1;
    if (e == 1) return a;

    if (logm == 0) logm = MSB(m);
    if (prec == 0) prec = BarrettPrecompute(m, logm);

    uint64_t res = 1; 

    a = a % m;

    while (e > 0) {
        if (e % 2 == 1) res = ModMultBarrett(res, a, m, prec, logm);
        //cout << e % 2 << " " << e << " " << a << " " << res << " " << m <<  endl;
        e >>= 1;
        a = ModMultBarrett(a, a, m, prec, logm);
        
    }

    return res;
}

bool IsPrime(uint64_t m, size_t iters) {
    // если n == 2 или n == 3 - эти числа простые, возвращаем true
    if (m == 2 || m == 3)
        return true;
 
    // если n < 2 или n четное - возвращаем false
    if (m < 2 || m % 2 == 0)
        return false;

    std::random_device rd;
    std::default_random_engine en(rd());
    std::uniform_int_distribution<int> rnd(2, m-2);


    uint64_t t = m - 1;

    uint64_t s = 0;
 
    while (t % 2 == 0) {
        t >>= 1;
        s += 1;
    }

    uint64_t logm = MSB(m);
    uint64_t prec = BarrettPrecompute(m, logm);

    // повторить k раз
    for (int i = 0; i < iters; i++) {
        // выберем случайное целое число a в отрезке [2, m − 2]
        uint64_t a = rnd(en);
 
        // x ← a^t mod n, вычислим с помощью возведения в степень по модулю
        uint64_t x = ModExp(a, t, m, prec, logm);
 
        // если x == 1 или x == m − 1, то перейти на следующую итерацию цикла
        if (x == 1 || x == m - 1)
            continue;
 
        // повторить s − 1 раз
        for (int r = 1; r < s; r++) {
            // x ← x^2 mod n
            x = ModMultBarrett(x, x, m, prec, logm);

            // если x == 1, то вернуть "составное"
            if (x == 1)
                return false;
 
            // если x == n − 1, то перейти на следующую итерацию внешнего цикла
            if (x == m - 1)
                break;
        }
 
        if (x != m - 1)
            return false;
    }
    
    return true;
}

uint64_t FindFirstPrimeDown(size_t logm, size_t M) {
    uint64_t m = ((static_cast<uint64_t>(1) << logm) + 1) - M;
    while (!IsPrime(m)) { 
        m = m - M;
    }
    return m;
}

uint64_t FindFirstPrimeUp(size_t logm, size_t M) {
    uint64_t m = ((static_cast<uint64_t>(1) << logm) + 1);
    while (!IsPrime(m)) { 
        m = m + M;
    }
    return m;
}

uint64_t FindPrevPrime(uint64_t m, size_t M) {
    uint64_t m1 = m - M;
    while (!IsPrime(m1)) { 
        m1 = m1 - M;
    }
    return m1;
}

uint64_t FindNextPrime(uint64_t m, size_t M) {
    uint64_t m1 = m + M;
    while (!IsPrime(m1)) { 
        m1 = m1 + M;
    }
    return m1;
}