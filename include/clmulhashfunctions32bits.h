/*
 * These are essentially scalar products functions.
 */
#ifndef CLMULHASFUNCTIONS32BITS_H_
#define CLMULHASFUNCTIONS32BITS_H_
#ifdef __PCLMUL__


#include "clmul.h"




// standard unoptimized 32-bit CLMUL hashing
uint32_t hashGaloisFieldMultilinear(const void *  randomsource, const uint32_t *  string, const size_t length) {
    const uint32_t * const endstring = string + length;
    const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    assert(sizeof(int) == sizeof(uint32_t));
    assert(((uintptr_t) randomsource32 & 15) == 0); // we expect cache line alignment for the keys
    __m128i acc = _mm_cvtsi32_si128(*(int * )(randomsource32 + length));
    for(; string!= endstring; ++randomsource32,++string ) {
        __m128i temp = _mm_set_epi64x(*randomsource32,*string);
        __m128i clprod  = _mm_clmulepi64_si128( temp, temp, 0x10);
        acc = _mm_xor_si128 (clprod,acc);
    }
    return barrettWithoutPrecomputation32(acc);
}


// optimized 32-bit CLMUL hashing
uint32_t hashGaloisFieldMultilinearHalfMultiplications(const void*  randomsource, const uint32_t *  string, const size_t length) {
    const uint32_t * const endstring = string + length;
    const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    assert(((uintptr_t) randomsource32 & 15) == 0); // we expect cache line alignment for the keys
    assert(sizeof(int) == sizeof(uint32_t));
    __m128i acc = _mm_cvtsi32_si128(*(int * )(randomsource32 + length));
    for(; string +3 < endstring; randomsource32+=4,string+=4 ) {
        const __m128i temp1 = _mm_load_si128((__m128i * )randomsource32);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i twosums = _mm_xor_si128(temp1,temp2);
        const __m128i part1 = _mm_unpacklo_epi32(twosums,_mm_setzero_si128());
        const __m128i clprod1  = _mm_clmulepi64_si128( part1, part1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        const __m128i part2 = _mm_unpackhi_epi32(twosums,_mm_setzero_si128());
        const __m128i clprod2  = _mm_clmulepi64_si128( part2, part2, 0x10);
        acc = _mm_xor_si128 (clprod2,acc);
    }
    if(string + 1  < endstring) {
        __m128i temp1 = _mm_set_epi64x(*randomsource32,*(randomsource32+1));
        __m128i temp2 = _mm_set_epi64x(*string,*(string+1));
        __m128i twosums = _mm_xor_si128(temp1,temp2);
        __m128i clprod  = _mm_clmulepi64_si128( twosums, twosums, 0x10);
        acc = _mm_xor_si128 (clprod,acc);
        randomsource32+=2;
        string+=2;
    }
    if( string!= endstring ) {
        __m128i temp = _mm_set_epi64x(*randomsource32,*string);
        __m128i clprod  = _mm_clmulepi64_si128( temp, temp, 0x10);
        acc = _mm_xor_si128 (clprod,acc);
    }
    return barrettWithoutPrecomputation32(acc);
}


#endif

#endif /* CLMULHASFUNCTIONS32BITS_H_ */
