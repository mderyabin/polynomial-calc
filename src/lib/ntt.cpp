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
}

NTT::~NTT() {
    if (tf) delete [] tf;
    if (itf) delete [] itf;
}


void NTT::ComputeForward(uint64_t *ax) {
    CooleyTukeyForwardNTT(ax, tf, N, m, prec, logm);
}

void NTT::ComputeForward(uint64_t *res, const uint64_t *ax){
    copy(ax, ax + N, res);
    CooleyTukeyForwardNTT(res, tf, N, m, prec, logm); // change function to avoid copy and gain performance
}

void NTT::ComputeInverse(uint64_t *ax){
    GentlemanSandeInverseNTT(ax, itf, N, m, invN, prec, prec_invN, logm);
}

void NTT::ComputeInverse(uint64_t *res, const uint64_t *ax) {
    copy(ax, ax + N, res);
    GentlemanSandeInverseNTT(res, itf, N, m, invN, prec, prec_invN, logm); // change function to avoid copy and gain performance
}

map<pair<size_t, uint64_t>, shared_ptr<NTT>> NTTManager::ntt_map = map<pair<size_t, uint64_t>, shared_ptr<NTT>>();
shared_ptr<NTT> NTTManager::GetNTTPtr(size_t N, uint64_t m) {
    if (N != (1<<(MSB(N)-1))) throw invalid_argument("dimension is invalid");
    if (!IsPrime(m) || (m % (2*N) != 1)) throw invalid_argument("modulus is invalid");

    pair<size_t, uint64_t> Nmpair(N, m);

    auto ntt_search = NTTManager::ntt_map.find(Nmpair);

    shared_ptr<NTT> ntt;
    if (ntt_search != NTTManager::ntt_map.end()) {
        ntt = ntt_search->second;
    }  else {
        ntt = make_shared<NTT>(N, m);
        pair<pair<size_t, uint64_t>, shared_ptr<NTT>> p(Nmpair, ntt);
        NTTManager::ntt_map.insert(p);
    }

    return ntt;
}

}