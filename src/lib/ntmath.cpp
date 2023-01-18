#include "ntmath.h"

#include <iostream>
#include <random>
using namespace std;

typedef unsigned __int128 uint128_t;

size_t MSB(uint64_t t) {
    size_t msb = 0;
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

uint64_t ModMult(uint64_t a, uint64_t b, uint64_t m) {
    size_t logm = MSB(m);
    uint64_t prec = BarrettPrecompute(m, logm);
    return ModMultBarrett(a, b, m, prec, logm);
}

void ModMultEq(uint64_t &a, uint64_t b, uint64_t m) {
    size_t logm = MSB(m);
    uint64_t prec = BarrettPrecompute(m, logm);
    ModMultBarrettEq(a, b, m, prec, logm);
}

//*
// https://core.ac.uk/download/pdf/287482281.pdf
uint64_t BarrettPrecompute(uint64_t m, size_t logm) {
    uint64_t bar = 2 * logm + 3;
    return (static_cast<uint128_t>(1) << bar) / m;
}

uint64_t ModMultBarrett(uint64_t a, uint64_t b, uint64_t m, uint64_t prec, size_t logm) {
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


void ModMultBarrettEq(uint64_t &a, uint64_t b, uint64_t m, uint64_t prec, size_t logm) {

    if (logm == 0) logm = MSB(m);
    if (prec == 0) prec = BarrettPrecompute(m, logm);

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


uint64_t ModExp(uint64_t a, uint64_t e, uint64_t m, uint64_t prec, size_t logm) {
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

    size_t logm = MSB(m);
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

uint64_t GCD(uint64_t a, uint64_t b) {
    while (a != b) {
        if (a > b) {
            uint64_t tmp = a;
            a = b;
            b = tmp;
        }
        b = b - a;
    }
    return a;
}

uint64_t abssub(uint64_t x, uint64_t y) {
    return (x > y) ? x - y : y - x;
}

uint64_t RhoPollard(uint64_t a) { 
    std::random_device rd;
    std::default_random_engine en(rd());
    std::uniform_int_distribution<int> rnd(2, a-2);

    uint64_t x = rnd(en);
    uint64_t y = 1;
    uint64_t i = 0;
    uint64_t stage = 2;

    size_t loga = MSB(a);
    uint64_t prec = BarrettPrecompute(a, loga);

    while (GCD(a, abssub(x, y)) == 1) {
        if (i == stage) {
            y = x;
            stage = stage * 2;
        }
        x = ModMultBarrett(x, x, a, prec, loga);
        x = ModAdd(x, 1, a);
        i = i + 1;
    } 
    return GCD(a, abssub(x, y));
}

vector<uint64_t> firstprimes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 
101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 
239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 
397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 
569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 
733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 
911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 
1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 
1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 
1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 
1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 
1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 
1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999};

vector<uint64_t> Factorize(uint64_t a) {
    vector<uint64_t> res;

    for (size_t i = 0; i < firstprimes.size(); i++) {
        uint64_t b = firstprimes[i];
        bool flag = true;
        while (flag) {
            if ( (a / b) * b == a ) {
                res.push_back(b);
                a /= b;
            } else {
                flag = false;
            }
        }
        if (a == 1) return res;
    }

    while (!IsPrime(a) && a != 1) {
        uint64_t b = RhoPollard(a);
        if (IsPrime(b)) 
            res.push_back(b);
        else {
            vector<uint64_t> bfact = Factorize(b);
            res.insert(res.end(), bfact.begin(), bfact.end());
        }
        a /= b;
    }

    if (a!=1) res.push_back(a);

    sort(res.begin(), res.end());

    return res;   
}

unordered_set<uint64_t> FindPrimefactors(uint64_t m) {
    unordered_set<uint64_t> res;
    vector<uint64_t> fact = Factorize(m);
    for (size_t i = 0; i < fact.size(); i++)
        res.insert(fact[i]);
    return res;
}

// Function to find smallest primitive root of n
uint64_t FindPrimitive(uint64_t n)
{
    unordered_set<uint64_t> s;
 
    // Check if n is prime or not
    if (!IsPrime(n))
        return 0;
 
    // Find value of Euler Totient function of n
    // Since n is a prime number, the value of Euler
    // Totient function is n-1 as there are n-1
    // relatively prime numbers.
    uint64_t phi = n-1;
 
    // Find prime factors of phi and store in a set
    s = FindPrimefactors(phi);

    size_t logn = MSB(n);
    uint64_t prec = BarrettPrecompute(n, logn);
 
    // Check for every number from 2 to phi
    for (int r=2; r<=phi; r++) {
        // Iterate through all prime factors of phi.
        // and check if we found a power with value 1
        bool flag = false;
        for (auto it = s.begin(); it != s.end(); it++) {
 
            // Check if r^((phi)/primefactors) mod n
            // is 1 or not
            if (ModExp(r, phi/(*it), n, prec, logn) == 1) {
                flag = true;
                break;
            }
         }
 
         // If there was no power with value 1.
         if (flag == false)
           return r;
    }
 
    // If no primitive root found
    return 0;
}

uint64_t FindGenerator(uint64_t m, size_t M) {
    uint64_t g0 = FindPrimitive(m);
    uint64_t g = ModExp(g0, (m - 1)/M, m);
    return g;
}

uint64_t ModInvPrime(uint64_t a, uint64_t m) {
    return ModExp(a, m-2, m);
}