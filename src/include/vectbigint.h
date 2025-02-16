#include <cmath>
#include <algorithm>

#include "ntt.h"

namespace polycalc {

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

uint64_t MulLazy(uint64_t *cx, const uint64_t *ax, uint64_t b, size_t n);
uint64_t MulLazy(uint64_t *cx, uint64_t b, size_t n);

void AddMod(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s);
void AddMod(uint64_t *cx, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s);

void MulMod(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s);
void MulMod(uint64_t *cx, const uint64_t *bx, const uint64_t *mx, size_t n, size_t s);

void SqrMod(uint64_t *cx, const uint64_t *ax, const uint64_t *mx, size_t n, size_t s);
void SqrMod(uint64_t *cx, const uint64_t *mx, size_t n, size_t s);

void FastModExp(uint64_t *cx, const uint64_t *ax, const uint64_t *ex, const uint64_t *mx, size_t n, size_t s);

bool IsBiggerOrEqual(const uint64_t *ax, const uint64_t *bx, size_t n);

void QMul(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, NTTInstance ntt);

// memory pool for convenience 
// class MemoryPool {
// public:
//     MemoryPool() {}
//     ~MemoryPool();

//     void Init(size_t n);

//     uint64_t *base;
//     uint64_t *shifted;
//     uint64_t *tmp1;
//     uint64_t *tmp2;
//     uint64_t *tmp3;

// };

// static MemoryPool globalMemPool;

}