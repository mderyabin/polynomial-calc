#include "encode.h"

using namespace std;

namespace polycalc {

const cdouble ZERO(0.0, 0.0);
const cdouble  ONE(1.0, 0.0);

cvector polyadd(cvector ax, cvector bx) { 
    size_t na  = ax.size();
    size_t nb  = bx.size();
    size_t nmax = (na > nb) ? na : nb;
    
    cvector cx(nmax);

    for (size_t i = 0; i < nmax; i++) {
        cdouble a = (i < na) ? ax[i] : ZERO;
        cdouble b = (i < nb) ? bx[i] : ZERO;
        cx[i] = a + b;
    }
    
    return cx;
}

cvector polymult(cvector ax, cvector bx) {
    size_t na  = ax.size();
    size_t nb  = bx.size();

    cvector cx(na + nb);

    for (size_t ij = 0; ij < na + nb; ij++) {
        cx[ij] = ZERO;
    }

    for (size_t i = 0; i < na; i++) {
        for (size_t j = 0; j < nb; j++) {
            cx[i + j] += ax[i] * bx[j];
        }
    }

    return cx;
}

// cpoly(x)
cdouble polyval(cvector cpoly, cdouble x) {
    cdouble xtmp = ONE;

    cdouble result = ZERO;

    for (size_t i = 0; i < cpoly.size(); i++) {
        result += cpoly[i] * xtmp;
        xtmp *= x;
    }
    
    return result;
}

cvector ComplexRootsOfUnity(size_t N) {
    cvector domain(N);

    for (size_t i = 0; i < N/2; i++){ 
        domain[i].real(cos(-2.0 * (2*i+1) * M_PI / (2.0*N)));
        domain[i].imag(sin(-2.0 * (2*i+1) * M_PI / (2.0*N)));

        domain[N - i - 1].real(cos(-2.0 * (2*i+1) * M_PI / (2.0*N)));
        domain[N - i - 1].imag(-sin(-2.0 * (2*i+1) * M_PI / (2.0*N)));
    }

    return domain;
}

// example can be found in https://eprint.iacr.org/2016/421.pdf (section 3.2)
Polynomial encode(cvector values, size_t N, uint64_t m, uint64_t scale) {
    size_t pow = values.size() * 2;

    Polynomial result(pow, m, true); 

    cvector domain = ComplexRootsOfUnity(N);

    cvector vals(pow);
    for (size_t i = 0; i < pow/2; i++) {
        vals[i] = values[i];
        vals[N - i - 1].real(values[i].real());
        vals[N - i - 1].imag(-values[i].imag());
    }

    cvector lagrange;
    cvector tmp, mono(2);

    tmp.reserve(pow);
    lagrange.reserve(pow);
    lagrange.resize(1);
    lagrange[0] = ZERO;
    cdouble c;
    for (size_t i = 0; i < pow; i++) {
        tmp.resize(1);
        tmp[0] = ONE;
        for (size_t j = 0; j < pow; j++) {
            if (i != j) {
                c = domain[i] - domain[j];
                mono[1] = ONE / c;
                mono[0] = - domain[j] / c;
                tmp = polymult(tmp, mono);
            }
        }
        for (size_t j = 0; j < tmp.size(); j++) {
            tmp[j] *= vals[i];
        }
        
        lagrange = polyadd(lagrange, tmp);
    }

    long double s = static_cast<long double>(scale);
    for (size_t i = 0; i < pow; i++) {
        long double a = round(s * lagrange[i].real());
        uint64_t ua = (a >= 0) ? static_cast<uint64_t>(a)%m : m - (static_cast<uint64_t>(-a))%m;
        result[i] = ua;
    }

    return result;
}

cvector decode(Polynomial poly, uint64_t scale) {
    size_t N = poly.GetN();
    uint64_t m = poly.GetModulus();

    long double dm = static_cast<long double>(m);
    long double ds = static_cast<long double>(scale);

    cvector result(N/2);

    cvector domain = ComplexRootsOfUnity(N);

    cvector cpoly(N);

    for (size_t i = 0; i < N; i++) {
        long double a = static_cast<long double>(poly[i]);
        a = (a < (dm/2)) ? a : (a - dm);
        a /= ds;
        cpoly[i].real(a);
    }
    
    for (size_t i = 0; i < N/2; i++) {
        result[i] = polyval(cpoly, domain[i]);
    }
    
    return result;
}

}