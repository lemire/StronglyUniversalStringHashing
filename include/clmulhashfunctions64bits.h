/*
 * clmulhashfunctions64bits.h
 *
 *  Created on: Jun 20, 2014
 *      Author: lemire
 */

#ifndef CLMULHASHFUNCTIONS64BITS_H_
#define CLMULHASHFUNCTIONS64BITS_H_
#ifdef __PCLMUL__

#include "clmul.h"



// a 64-bit version
// strongly universal and regular
uint64_t hashGaloisFieldfast64(const void*  rs, const uint64_t *  string, const size_t length) {
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )rs;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string+1< endstring; randomsource+=2,string+=2 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        const __m128i clprod2  = _mm_clmulepi64_si128( temp1, temp2, 0x11);
        acc = _mm_xor_si128 (clprod1,acc);
        acc = _mm_xor_si128 (clprod2,acc);
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_set_epi64x(0,*randomsource);
        const __m128i temp2 = _mm_set_epi64x(0,*string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return barrettWithoutPrecomputation64(acc);
}

uint64_t hashGaloisFieldfast64_precomp(const void*  rs, const uint64_t *  string, const size_t length) {
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )rs;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string+1< endstring; randomsource+=2,string+=2 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        const __m128i clprod2  = _mm_clmulepi64_si128( temp1, temp2, 0x11);
        acc = _mm_xor_si128 (clprod1,acc);
        acc = _mm_xor_si128 (clprod2,acc);
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_set_epi64x(0,*randomsource);
        const __m128i temp2 = _mm_set_epi64x(0,*string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return precompReduction64(acc);
}

// a 64-bit version
// strongly universal and regular
uint64_t hashGaloisFieldfast64_precomp_unroll(const void*  rs, const uint64_t *  string, const size_t length) {
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )rs;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string+1< endstring; randomsource+=4,string+=4 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        const __m128i clprod2  = _mm_clmulepi64_si128( temp1, temp2, 0x11);
        acc = _mm_xor_si128 (clprod1,acc);
        acc = _mm_xor_si128 (clprod2,acc);
        const __m128i temp1_2 = _mm_lddqu_si128((__m128i * )(randomsource + 2));
        const __m128i temp2_2 = _mm_lddqu_si128((__m128i *) (string + 2));
        const __m128i clprod1_2  = _mm_clmulepi64_si128( temp1_2, temp2_2, 0x00);
        const __m128i clprod2_2  = _mm_clmulepi64_si128( temp1_2, temp2_2, 0x11);
        acc = _mm_xor_si128 (clprod1_2,acc);
        acc = _mm_xor_si128 (clprod2_2,acc);

    }
    for(; string+1< endstring; randomsource+=2,string+=2 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        const __m128i clprod2  = _mm_clmulepi64_si128( temp1, temp2, 0x11);
        acc = _mm_xor_si128 (clprod1,acc);
        acc = _mm_xor_si128 (clprod2,acc);
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_set_epi64x(0,*randomsource);
        const __m128i temp2 = _mm_set_epi64x(0,*string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return precompReduction64(acc);
}


// a 64-bit version with half the number of multiplications
// strongly universal but not regular
uint64_t hashGaloisFieldfast64half(const void*  rs, const uint64_t *  string, const size_t length) {
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )rs;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string+1< endstring; randomsource+=2,string+=2 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_set_epi64x(0,*randomsource);
        const __m128i temp2 = _mm_set_epi64x(0,*string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return barrettWithoutPrecomputation64(acc);
}

// a 64-bit version with half the number of multiplications
// strongly universal but not regular
uint64_t hashGaloisFieldfast64half_precomp(const void*  rs, const uint64_t *  string, const size_t length) {
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )rs;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string+1< endstring; randomsource+=2,string+=2 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_set_epi64x(0,*randomsource);
        const __m128i temp2 = _mm_set_epi64x(0,*string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return precompReduction64(acc);
}


// a 64-bit version with half the number of multiplications
uint64_t hashGaloisFieldfast64halfunrolled(const void*  rs, const uint64_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )rs;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string+3 < endstring; randomsource+=4,string+=4 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        const __m128i temp12 = _mm_lddqu_si128((__m128i * )(randomsource + 2));
        const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string+2));
        const __m128i add12 =  _mm_xor_si128 (temp12,temp22);
        const __m128i clprod12  = _mm_clmulepi64_si128( add12, add12, 0x10);
        acc = _mm_xor_si128 (clprod12,acc);
    }
    if(string+1< endstring) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        randomsource+=2;
        string+=2;
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_set_epi64x(0,*randomsource);
        const __m128i temp2 = _mm_set_epi64x(0,*string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return barrettWithoutPrecomputation64(acc);
}

// a 64-bit version with half the number of multiplications
uint64_t hashGaloisFieldfast64halfunrolled_precomp(const void*  rs, const uint64_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )rs;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string+3 < endstring; randomsource+=4,string+=4 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        const __m128i temp12 = _mm_lddqu_si128((__m128i * )(randomsource + 2));
        const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string+2));
        const __m128i add12 =  _mm_xor_si128 (temp12,temp22);
        const __m128i clprod12  = _mm_clmulepi64_si128( add12, add12, 0x10);
        acc = _mm_xor_si128 (clprod12,acc);
    }
    if(string+1< endstring) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        randomsource+=2;
        string+=2;
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_set_epi64x(0,*randomsource);
        const __m128i temp2 = _mm_set_epi64x(0,*string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return precompReduction64(acc);
}

#endif

#endif /* CLMULHASHFUNCTIONS64BITS_H_ */
