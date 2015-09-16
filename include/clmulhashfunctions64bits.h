/*
 * These are essentially scalar products functions.
 */

#ifndef CLMULHASHFUNCTIONS64BITS_H_
#define CLMULHASHFUNCTIONS64BITS_H_
#ifdef __PCLMUL__

#include "clmul.h"


/// These are multilinear functions. Fast but maybe not practical. Just used for comparison purposes.

// a 64-bit multilinear  function
// strongly universal
uint64_t hashGaloisFieldfast64_precomp_unroll(const void* rs, const uint64_t * string, const size_t length) {
    const uint64_t * const endstring = string + length;
    const uint64_t * randomsource = ( const uint64_t * )rs;
    assert(((uintptr_t) randomsource & 15) == 0); // we expect cache line alignment for the keys
    __m128i acc = _mm_loadl_epi64((__m128i const*)(randomsource + length));
    for(; string+3 < endstring; randomsource+=4,string+=4 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_load_si128((__m128i *) string);
        const __m128i clprod1 = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        const __m128i clprod2 = _mm_clmulepi64_si128( temp1, temp2, 0x11);
        acc = _mm_xor_si128 (clprod1,acc);
        acc = _mm_xor_si128 (clprod2,acc);
        const __m128i temp12 = _mm_load_si128((__m128i * )(randomsource + 2));
        const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string+2));
        const __m128i clprod12 = _mm_clmulepi64_si128( temp12, temp22, 0x00);
        const __m128i clprod22 = _mm_clmulepi64_si128( temp12, temp22, 0x11);
        acc = _mm_xor_si128 (clprod12,acc);
        acc = _mm_xor_si128 (clprod22,acc);
    }
    if(string+1< endstring) {
        const __m128i temp1 = _mm_load_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i clprod1 = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        const __m128i clprod2 = _mm_clmulepi64_si128( temp1, temp2, 0x11);
        acc = _mm_xor_si128 (clprod1,acc);
        acc = _mm_xor_si128 (clprod2,acc);
        randomsource+=2;
        string+=2;
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_load_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_loadl_epi64((__m128i const*)string);
        const __m128i clprod1 = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return precompReduction64(acc);
}

// a 64-bit multilinear with half the number of multiplications
uint64_t hashGaloisFieldfast64halfunrolled_precomp(const void* rs, const uint64_t * string, const size_t length) {
    const uint64_t * const endstring = string + length;
    const uint64_t * randomsource = ( const uint64_t * )rs;
    assert(((uintptr_t) randomsource & 15) == 0); // we expect cache line alignment for the keys
    __m128i acc = _mm_loadl_epi64((__m128i const*)(randomsource + length));
    for(; string+3 < endstring; randomsource+=4,string+=4 ) {
        const __m128i temp1 = _mm_load_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 = _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        const __m128i temp12 = _mm_load_si128((__m128i * )(randomsource + 2));
        const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string+2));
        const __m128i add12 = _mm_xor_si128 (temp12,temp22);
        const __m128i clprod12 = _mm_clmulepi64_si128( add12, add12, 0x10);
        acc = _mm_xor_si128 (clprod12,acc);
    }
    if(string+1< endstring) {
        const __m128i temp1 = _mm_load_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 = _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        randomsource+=2;
        string+=2;
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_load_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_loadl_epi64((__m128i const*)string);
        const __m128i clprod1 = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return precompReduction64(acc);
}

#endif

#endif /* CLMULHASHFUNCTIONS64BITS_H_ */
