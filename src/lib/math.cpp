#include "../include/math.h"

uint64_t ModAdd(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t c = a + b;
    return (c >= m) ? c - m : c;
}