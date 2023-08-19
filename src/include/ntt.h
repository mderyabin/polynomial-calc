/**
 * @file ntt.h
 * @author Maxim Deryabin (maxim.deryabin@gmail.com)
 * @brief Organization of the NTT/INTT.
 * @version 0.1
 * @date 2023-05-14
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef __NTT_H__
#define __NTT_H__

#include <memory>
#include <map>

#include "polymath.h"

namespace polycalc {

/**
 * @brief High level wrapper for NTT/INTT transformations.
 * 
 * Allows to organize the computations for fixed prime modulus. 
 * 
 */
class NTT {
private:
    size_t N;           // ring dimensions
    size_t logN;
    uint64_t m;         // modulus

    uint64_t invN;      // inversion of N

    uint64_t *tf;       // twiddle factors for forward transform
    uint64_t *itf;      // twiddle factors for inverse transform

    uint64_t *prec_tf;  // constants precomputed for Shoup mult for each tf
    uint64_t *prec_itf; // constants precomputed for Shoup mult for each itf

    size_t logm;        // bitsize of modulus
    uint64_t prec;      // precomputed constant for Barrett multiplication
    uint64_t prec_invN; // precomputed constant for Shoup multiplication

public:
    NTT(size_t _N, uint64_t _m);
    ~NTT();

    /**
     * @brief In-place forward NTT.
     * 
     * @param ax input polynomial as native C array (should be in COEF format)
     */
    void ComputeForward(uint64_t *ax);

    /**
     * @brief Forward NTT.
     * 
     * @param res output polynomial as native C array
     * @param ax input polynomial as native C array (should be in COEF format)
     */
    void ComputeForward(uint64_t *res, const uint64_t *ax);

    /**
     * @brief In-place inverse NTT.
     * 
     * @param ax input polynomial as native C array (should be in EVAL format)
     */
    void ComputeInverse(uint64_t *ax);

    /**
     * @brief Inverse NTT.
     * 
     * @param res output polynomial as native C array
     * @param ax input polynomial as native C array (should be in EVAL format)
     */
    void ComputeInverse(uint64_t *res, const uint64_t *ax);
};


/**
 * @brief Manager for NTT classes for different moduli.
 * 
 * Provides dynamic organization for existing NTT classes for each modulus to avoid re-declaration
 * and re-initialization of NTT class for same modulus.
 */
class NTTManager {
private:

    // Internal map of NTT classes.
    // Each entity in the map corresponds to particular pair N and m.
    static std::map<std::pair<size_t, uint64_t>, std::shared_ptr<NTT>> ntt_map;
public:

    /**
     * @brief Managing the map of NTT classes. 
     * 
     * Preserve singular instance of NTT for each pair N and m. 
     * Allows to avoid re-declaration of the NTT object with same parameters.
     * 
     * @param N dimension
     * @param m modulus
     * @return shared pointer to the objet of NTT class
     */
    static std::shared_ptr<NTT> GetNTTPtr(size_t N, uint64_t m);
};

}

#endif /* __NTT_H__ */