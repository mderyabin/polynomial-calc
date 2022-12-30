#include "polynomial.h"

#include <algorithm>

using namespace std;

Polynomial::Polynomial(size_t _N, uint64_t _m, bool initWithZeros) : N(_N), m(_m) {
    ax = new uint64_t[N];
    if (initWithZeros) {
        for (size_t i = 0; i < N; i++) {
            ax[i] = 0L;
        }
    }
    logm = MSB(m);
    mu = BarrettPrecompute(m, logm);
}

Polynomial::Polynomial(uint64_t *_ax, size_t _N, uint64_t _m) : N(_N), m(_m) {
    ax = new uint64_t[N];
    copy(_ax, _ax + N, ax);
    logm = MSB(m);
    mu = BarrettPrecompute(m, logm);
}

Polynomial::Polynomial(const Polynomial &o) : N(o.N), m(o.m), logm(o.logm), mu(o.mu) {
    ax = new uint64_t[N];
    copy(o.ax, o.ax + N, ax);
}

Polynomial& Polynomial::operator=(const Polynomial &o) {
    N = o.N;
    m = o.m;
    logm = o.logm;
    mu = o.mu;
    ax = new uint64_t[N];
    copy(o.ax, o.ax + N, ax);

    return *this;
}

Polynomial::Polynomial(Polynomial &&o) : N(o.N), m(o.m), logm(o.logm), mu(o.mu) {
    move(o.ax, o.ax + N, ax);

    o.N = 0;
    o.m = 0;
    o.logm = 0;
    o.mu = 0;
    o.ax = NULL;
}

Polynomial& Polynomial::operator=(Polynomial &&o) {
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

uint64_t Polynomial::operator()(uint64_t x) const {
    uint64_t res = ax[0];
    uint64_t xpow = x;
    uint64_t pr = ShoupPrecompute(x, m);
    for (size_t i = 1; i < N; i++) {
        res = ModAdd(res, ModMultBarrett(ax[i], xpow, m, mu, logm), m);
        xpow = ModMulShoup(xpow, x, m, logm);
    }
    
    return res;
}

ostream& operator<<(ostream& os, const Polynomial& Polynomial) {
    size_t N = Polynomial.GetN();
    
    if (N >=1) os << Polynomial.ax[0];

    for (size_t i = 1; i < N; i++) {
        os << " + " << Polynomial.ax[i] << "X^" << i;
    }
    
    return os;
}

const Polynomial& operator+=(Polynomial& left, const Polynomial& right) {
    if (left.m != right.m || left.m != right.m)
        return left;
    ModAdd(left.ax, right.ax, left.N, left.m);
    return left;
}

const Polynomial operator+(const Polynomial& left, const Polynomial& right) {
    Polynomial res(left);
    return (res += right);
}

const Polynomial operator*(const Polynomial& left, const Polynomial& right) {
    Polynomial res(left);
    NaiveNegacyclicConvolution(res.ax, left.ax, right.ax, res.N, res.m, res.mu, res.logm);
    return res;
}

const Polynomial& operator*=(Polynomial& left, const Polynomial& right) {
    uint64_t *temp = new uint64_t[left.N];
    NaiveNegacyclicConvolution(temp, left.ax, right.ax, left.N, left.m, left.mu, left.logm);
    move(temp, temp + left.N, left.ax);
    delete [] temp;
    return left;
}


