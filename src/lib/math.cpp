#include "../include/math.h"

#include <iostream>
using namespace std;

uint64_t ModAdd(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t c = a + b;
    return (c >= m) ? c - m : c;
}

//*
// https://core.ac.uk/download/pdf/287482281.pdf
uint64_t BarrettPrecompute(uint64_t m, uint64_t logm) {
    uint64_t bar = 2 * logm;
    return (static_cast<unsigned __int128>(1) << bar) / m;
}

uint64_t ModMultBarrett(uint64_t a, uint64_t b, uint64_t m, uint64_t prec, uint64_t logm) {
	uint64_t res;
    //int64_t bar = 2 * logm;

    uint64_t mul = a * b;
    uint64_t tmp = mul >> (logm-1);
    tmp *= prec;
    tmp >>= (logm + 1);

    tmp *= m;
    res = (mul - tmp) & ( (1 << (logm + 1)) - 1 ); // mod 2 ^(logm + 1)
    
    while (res >= m) res -= m;

    return res;
}
//*/

/* то же самое, только 64-бит версия из HEAAN 
// https://github.com/snucrypto/HEAAN/blob/master/HEAAN/src/RingMultiplier.cpp

// logm > 32 !!
uint64_t BarrettPrecompute(uint64_t m, uint64_t logm) {
    uint64_t bar = 2 * logm;
    return (static_cast<unsigned __int128>(1) << bar) / m;
}

// logm > 32 !!
uint64_t ModMultBarrett(uint64_t a, uint64_t b, uint64_t m, uint64_t prec, uint64_t logm) {
	uint64_t res;
    uint64_t bar = 2 * logm;
    
    unsigned __int128 mul = static_cast<unsigned __int128>(a) * b;
	uint64_t abot = static_cast<uint64_t>(mul);
	uint64_t atop = static_cast<uint64_t>(mul >> 64);
	unsigned __int128 tmp = static_cast<unsigned __int128>(abot) * prec;
	tmp >>= 64;
	tmp += static_cast<unsigned __int128>(atop) * prec;
	tmp >>= bar - 64;
	tmp *= m;
	tmp = mul - tmp;
	res = static_cast<uint64_t>(tmp);  // mod 2^64
	if (res >= m) res -= m;

    return res;
}

//*/