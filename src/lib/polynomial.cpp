#include "polynomial.h"

#include <algorithm>


using namespace std;

namespace polycalc {

Polynomial::Polynomial(size_t _N, uint64_t _m, bool initWithZeros, Format _format) : N(_N), m(_m), format(_format) {
    if (N != (1<<(MSB(N)-1))) throw invalid_argument("dimension is invalid");
    if (!IsPrime(m) || (m % (2*N) != 1)) throw invalid_argument("modulus is invalid");
    
    ax = new uint64_t[N];
    if (initWithZeros) 
        SetZero();

    logm = MSB(m);
    mu = BarrettPrecompute(m, logm);

    ntt = NTTManager::GetNTTPtr(N, m);
}

Polynomial::Polynomial(uint64_t *_ax, size_t _N, uint64_t _m, Format _format) : N(_N), m(_m), format(_format) {
    if (N != (1<<(MSB(N)-1))) throw invalid_argument("dimension is invalid");
    if (!IsPrime(m) && (m % (2*N) != 1)) throw invalid_argument("modulus is invalid");

    ax = new uint64_t[N];
    copy(_ax, _ax + N, ax);
    logm = MSB(m);
    mu = BarrettPrecompute(m, logm);

    ntt = NTTManager::GetNTTPtr(N, m);
}

Polynomial::Polynomial(const Polynomial &o) : N(o.N), m(o.m), logm(o.logm), mu(o.mu), format(o.format) {
    ax = new uint64_t[N];
    ntt = o.ntt;
    copy(o.ax, o.ax + N, ax);
}

Polynomial& Polynomial::operator=(const Polynomial &o) {
    N = o.N;
    m = o.m;
    logm = o.logm;
    mu = o.mu;
    format = o.format;
    ax = new uint64_t[N];
    ntt = o.ntt;
    copy(o.ax, o.ax + N, ax);

    return *this;
}

Polynomial::Polynomial(Polynomial &&o) : N(o.N), m(o.m), logm(o.logm), mu(o.mu), format(o.format) {
    ntt = o.ntt;
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
    format = o.format;
    ntt = o.ntt;    // do not delete ntt! 

    o.N = 0;
    o.m = 0;
    o.logm = 0;
    o.mu = 0;
    o.ax = NULL;

    return *this;
}

void Polynomial::GenerateUniform(Format _format) {
    format = _format;
    GenerateUniform();
}

void Polynomial::GenerateUniform() {
    GenerateUniformPoly(ax, N, m);  
}

void Polynomial::SetZero(Format _format) {
    format = _format;
    SetZero();
}

void Polynomial::SetZero() {
    for (size_t i = 0; i < N; i++) {
        ax[i] = 0L;
    }
}

uint64_t Polynomial::operator()(uint64_t x) const {
    x %= m;
    if (format == COEF)
        return ComputeValue(x, ax, N, m);
    else {
        uint64_t *coeffs = new uint64_t[N]; 
        copy(ax, ax + N, coeffs);

        // todo: transform coeffs using INTT

        uint64_t res = ComputeValue(x, coeffs, N, m);

        delete [] coeffs;

        return res;
    } 
}

ostream& operator<<(ostream& os, const Polynomial& polynomial) {
    size_t N = polynomial.GetN();
	uint64_t *ax = polynomial.ax;
    
    if (polynomial.format == COEF) {
        bool printplus = false;
        if (N >=1 ) {
            for (size_t i = 0; i < N; i++) {
                if (ax[i] != 0) {
                    if (printplus) os << " + ";
                    else printplus = true;
                    if (ax[i] != 1 || i == 0) os << ax[i]; 
                    if (i != 0) os << "X^" << i;
                }
            }
        } else {
            os << 0;
        }
    } else {
        os << " [" << ax[0];
        for (size_t i = 1; i < N; i++) {
            os << ", " << ax[i];
        }
        os << "]";
    }

    return os;
}

void Polynomial::SetFormatEval() {
    if (format == COEF) {
        (*ntt).ComputeForward(ax);
        format = EVAL;
    }
}

void Polynomial::SetFormatCoef() {
    if (format == EVAL) {
        (*ntt).ComputeInverse(ax);
        format = COEF;
    }
}

const Polynomial& operator+=(Polynomial& left, const Polynomial& right) {
    if (left.m != right.m || left.N != right.N)
        throw runtime_error("Parameters mismatched in addition");
    ModAdd(left.ax, right.ax, left.N, left.m);
    return left;
}

const Polynomial operator+(const Polynomial& left, const Polynomial& right) {
    if (left.m != right.m || left.N != right.N)
        throw runtime_error("Parameters mismatched in addition");
    Polynomial res(left);
    return (res += right);
}

const Polynomial operator*(const Polynomial& left, const Polynomial& right) {
    if (left.m != right.m || left.N != right.N)
        throw runtime_error("Parameters mismatched in multiplication");

    Polynomial res(left);
    if (left.format == COEF && right.format == COEF)
        NaiveNegacyclicConvolution(res.ax, left.ax, right.ax, res.N, res.m, res.mu, res.logm);
    else if (left.format == EVAL && right.format == EVAL)
        ModHadamardMul(res.ax, right.ax, res.N, res.m, res.mu, res.logm);
    else 
        throw runtime_error("Formats mismatched in multiplication");
    return res;
}

const Polynomial& operator*=(Polynomial& left, const Polynomial& right) {
    if (left.m != right.m || left.N != right.N)
        throw runtime_error("Parameters mismatched in multiplication");

    if (left.format == COEF && right.format == COEF) {
        uint64_t *temp = new uint64_t[left.N];
        NaiveNegacyclicConvolution(temp, left.ax, right.ax, left.N, left.m, left.mu, left.logm);
        move(temp, temp + left.N, left.ax);
        delete [] temp;
    } else if (left.format == EVAL && right.format == EVAL) {
        ModHadamardMul(left.ax, right.ax, left.N, left.m, left.mu, left.logm);
    } else 
        throw runtime_error("Formats mismatched in multiplication");
    return left;
}

// order: format (0 for EVAL and 1 for COEF), N, m, logm, mu, ax[i]
template<class Archive>
void Polynomial::save(Archive & archive) const {

    // step 1: save format
    int format_to_save = 0;
    
    switch (format) {
    case EVAL: 
        format_to_save = 0;
        break;
    case COEF:
    default:
        format_to_save = 1;
        break;
    }

    archive(cereal::make_nvp("f", format_to_save));

    // step 2: save N, m, logm, mu
    archive(cereal::make_nvp("N", N));
    archive(cereal::make_nvp("m", m));
    archive(cereal::make_nvp("l", logm));
    archive(cereal::make_nvp("u", mu));

    // step 3: save ax
    for (size_t i = 0; i < N; i++) {
        archive(ax[i]);
    }
    
}

template<class Archive>
void Polynomial::load(Archive & archive) {

    // step 1: load format
    int format_to_load = 0;
    archive(cereal::make_nvp("f", format_to_load));

    switch (format_to_load) {
    case 0: 
        format = EVAL;
        break;
    case 1:
    default:
        format = COEF;
        break;
    }

    // step 2: load N, m, logm, mu
    size_t N_temp;
    archive(cereal::make_nvp("N", N_temp));
    archive(cereal::make_nvp("m", m));
    archive(cereal::make_nvp("l", logm));
    archive(cereal::make_nvp("u", mu));

    if (N_temp != N) {
        N = N_temp;
        if (ax) delete [] ax;
        ax = NULL;
    }

    if (!ax) ax = new uint64_t[N];

    // step 3: load ax
    for (size_t i = 0; i < N; i++) {
        archive(ax[i]);
    }
}

void Polynomial::Serialize(string filename, SER_Archive_Type TYPE) {
    SerializeObject(*this, filename, TYPE);
}

void Polynomial::Deserialize(std::string filename, SER_Archive_Type TYPE) {
    DeserializeObject(*this, filename, TYPE);
}

}