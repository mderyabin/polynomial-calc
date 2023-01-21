#ifndef __NTT_H__
#define __NTT_H__

#include <memory>
#include <map>

#include "polymath.h"


class NTT {
private:
    size_t N;           // ring dimensions
    uint64_t m;         // modulus

    uint64_t invN;      // inversion of N

    uint64_t *tf;       // twiddle factors for forward transform
    uint64_t *itf;      // twiddle factors for inverse transform

    size_t logm;        // bitsize of modulus
    uint64_t prec;      // precomputed constant for Barrett multiplication
    uint64_t prec_invN; // precomputed constant for Shoup multiplication

public:
    NTT(size_t _N, uint64_t _m);
    ~NTT();

    void ComputeForward(uint64_t *ax);
    void ComputeForward(uint64_t *res, const uint64_t *ax);

    void ComputeInverse(uint64_t *ax);
    void ComputeInverse(uint64_t *res, const uint64_t *ax);
};



class NTTManager {
    static std::map<std::pair<size_t, uint64_t>, std::shared_ptr<NTT>> ntt_map;
public:
    static std::shared_ptr<NTT> GetNTTPtr(size_t N, uint64_t m);
};

#endif /* __NTT_H__ */