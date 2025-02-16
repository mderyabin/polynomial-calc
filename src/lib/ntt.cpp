#include "ntt.h"

#include <iostream>

using namespace std;

namespace polycalc {

NTT::NTT(size_t _N, uint64_t _m) : N(_N), m(_m) {
    if (N != (1<<(MSB(N)-1))) throw invalid_argument("dimension is invalid");
    if (!IsPrime(m) || (m % (2*N) != 1)) throw invalid_argument("modulus is invalid");

    tf = new uint64_t[N];
    itf = new uint64_t[N];

    logm = MSB(m);
    logN = MSB(N)-1;

    invN = ModInvPrime(N, m);

    prec = BarrettPrecompute(m, logm);
    prec_invN = ShoupPrecompute(invN, m);

    // compute twiddle factors for forward NTT
    // combines powers of N-th root of unity for NTT with powers of 2N-th root for NWC
    // stored in bit-reversed order
    ComputeTwiddleFactors(tf, N, m, false);

    // compute twiddle factors for inverse NTT
    // combines powers of inversion of N-th root of unity for INTT with powers of inversion of 2N-th root for INWC
    // stored in bit-reversed order
    ComputeTwiddleFactors(itf, N, m, true);

    prec_tf = new uint64_t[N];
    prec_itf = new uint64_t[N];

    ShoupPrecompute(prec_tf, tf, N, m);
    ShoupPrecompute(prec_itf, itf, N, m);
}

NTT::~NTT() {
    if (tf) delete [] tf;
    if (itf) delete [] itf;
    if (itf) delete [] prec_tf;
    if (itf) delete [] prec_itf;
}


void NTT::ComputeForward(uint64_t *ax) {
    // CooleyTukeyForwardNTT(ax, tf, N, m, prec, logm);
    CooleyTukeyForwardNTT(ax, tf, N, m, prec_tf, logN);
}

void NTT::ComputeForward(uint64_t *res, const uint64_t *ax){
    copy(ax, ax + N, res);
    //CooleyTukeyForwardNTT(res, tf, N, m, prec, logm); // change function to avoid copy and gain performance
    CooleyTukeyForwardNTT(res, tf, N, m, prec_tf, logN);
}

void NTT::ComputeInverse(uint64_t *ax){
    // GentlemanSandeInverseNTT(ax, itf, N, m, invN, prec, prec_invN, logm);
    GentlemanSandeInverseNTT(ax, itf, N, m, invN, prec_itf, prec_invN);
}

void NTT::ComputeInverse(uint64_t *res, const uint64_t *ax) {
    copy(ax, ax + N, res);
    // GentlemanSandeInverseNTT(res, itf, N, m, invN, prec, prec_invN, logm); // change function to avoid copy and gain performance
    GentlemanSandeInverseNTT(res, itf, N, m, invN, prec_itf, prec_invN);
}

map<pair<size_t, uint64_t>, NTTInstance> NTTManager::ntt_map = map<pair<size_t, uint64_t>, NTTInstance>();
NTTInstance NTTManager::GetNTTPtr(size_t N, uint64_t m) {
    if (N != (1<<(MSB(N)-1))) throw invalid_argument("dimension is invalid");
    if (!IsPrime(m) || (m % (2*N) != 1)) throw invalid_argument("modulus is invalid");

    pair<size_t, uint64_t> Nmpair(N, m);

    auto ntt_search = NTTManager::ntt_map.find(Nmpair);

    NTTInstance ntt;
    if (ntt_search != NTTManager::ntt_map.end()) {
        ntt = ntt_search->second;
    }  else {
        ntt = make_shared<NTT>(N, m);
        pair<pair<size_t, uint64_t>, NTTInstance> p(Nmpair, ntt);
        NTTManager::ntt_map.insert(p);
    }

    return ntt;
}

}