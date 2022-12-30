#include "poly.h"

#include <algorithm>

using namespace std;

Poly::Poly(size_t _N, uint64_t _m, bool initWithZeros) : N(_N), m(_m) {
    ax = new uint64_t[N];
    if (initWithZeros) {
        for (size_t i = 0; i < N; i++) {
            ax[i] = 0L;
        }
    }
    logm = MSB(m);
    mu = BarrettPrecompute(m, logm);
}

Poly::Poly(uint64_t *_ax, size_t _N, uint64_t _m) : N(_N), m(_m) {
    ax = new uint64_t[N];
    copy(_ax, _ax + N, ax);
    logm = MSB(m);
    mu = BarrettPrecompute(m, logm);
}

Poly::Poly(const Poly &o) : N(o.N), m(o.m), logm(o.logm), mu(o.mu) {
    ax = new uint64_t[N];
    copy(o.ax, o.ax + N, ax);
}

Poly& Poly::operator=(const Poly &o) {
    N = o.N;
    m = o.m;
    logm = o.logm;
    mu = o.mu;
    ax = new uint64_t[N];
    copy(o.ax, o.ax + N, ax);

    return *this;
}

Poly::Poly(Poly &&o) : N(o.N), m(o.m), logm(o.logm), mu(o.mu) {
    move(o.ax, o.ax + N, ax);

    o.N = 0;
    o.m = 0;
    o.logm = 0;
    o.mu = 0;
    o.ax = NULL;
}

Poly& Poly::operator=(Poly &&o) {
    move(o.ax, o.ax + N, ax);
    N = o.N;
    m = o.m;
    logm = o.logm;
    mu = o.mu;

    o.N = 0;
    o.m = 0;
    o.logm = 0;
    o.mu = 0;
    o.ax = NULL;

    return *this;
}

uint64_t Poly::operator()(uint64_t x) const {
    uint64_t res = ax[0];
    uint64_t xpow = x;
    uint64_t pr = ShoupPrecompute(x, m);
    for (size_t i = 1; i < N; i++) {
        res = ModAdd(res, ModMultBarrett(ax[i], xpow, m, mu, logm), m);
        xpow = ModMulShoup(xpow, x, m, logm);
    }
    
    return res;
}

ostream& operator<<(ostream& os, const Poly& poly) {
    size_t N = poly.GetN();
    
    if (N >=1) os << poly.ax[0];

    for (size_t i = 1; i < N; i++) {
        os << " + " << poly.ax[i] << "X^" << i;
    }
    
    return os;
}

const Poly& operator+=(Poly& left, const Poly& right) {
    if (left.GetN() != right.GetN() || left.GetModulus() != right.GetModulus())
        return left;
    for (size_t i = 0; i < left.GetN(); i++) {
        left[i] = ModAdd(left[i], right.at(i), left.GetModulus());
    }
    return left;
}

const Poly operator+(const Poly& left, const Poly& right) {
    Poly res(left);
    return (res += right);
}
