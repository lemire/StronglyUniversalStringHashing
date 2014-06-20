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
uint64_t hashGaloisFieldfast64halfunrolled(const void*  rs, const uint64_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )randomsource;
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


// simple 64-bit polynomial hashing, uses only one key
// not expected to be fast!
uint64_t hashGaloisFieldPoly64(const void*  rs, const uint64_t *  string, const size_t length) {
    const uint64_t *  randomsource = ( const uint64_t * )randomsource;
	assert(*randomsource != 0);//otherwise silly
    const uint64_t * const endstring = string + length;
    __m128i key = _mm_set_epi64x(0,*(randomsource));
    __m128i acc = _mm_set_epi64x(0,*string);
    ++string;
    for(; string  < endstring; ++string ) {
        const __m128i temp = _mm_set_epi64x(0,*string);
        const __m128i multi = _mm_clmulepi64_si128( acc, key, 0x00);
        acc = barrettWithoutPrecomputation64_si128(multi);
        acc = _mm_xor_si128 (acc,temp);
    }
    __m128i multi = _mm_clmulepi64_si128( acc, key, 0x00);
    return barrettWithoutPrecomputation64(multi);
}

uint64_t precomphashGaloisFieldPoly64(const void*  rs, const uint64_t *  string, const size_t length) {
	const uint64_t *  randomsource = ( const uint64_t * )randomsource;
	assert(*randomsource != 0);//otherwise silly
    const uint64_t * const endstring = string + length;
    __m128i key = _mm_set_epi64x(0,*(randomsource));
    __m128i acc = _mm_set_epi64x(0,*string);
    ++string;
    for(; string  < endstring; ++string ) {
        const __m128i temp = _mm_set_epi64x(0,*string);
        const __m128i multi = _mm_clmulepi64_si128( acc, key, 0x00);
        acc = precompReduction64_si128(multi);
        acc = _mm_xor_si128 (acc,temp);
    }
    __m128i multi = _mm_clmulepi64_si128( acc, key, 0x00);
    return precompReduction64(multi);
}

// fast 64-bit polynomial hashing, uses only one key
// expected to be fast!
//TODO: can use more keys for increased universality
uint64_t fasthashGaloisFieldPoly64(const void*  rs, const uint64_t *  string, const size_t length) {
	const uint64_t *  randomsource = ( const uint64_t * )rs;
	assert(*randomsource != 0);//otherwise silly
    const uint64_t * const endstring = string + length;
    __m128i tkey1 = _mm_set_epi64x(0,*(randomsource));
    // we start by precomputing the powers of the key
    __m128i tkey2 = barrettWithoutPrecomputation64_si128(
       _mm_clmulepi64_si128( tkey1, tkey1, 0x00));
    __m128i tkey3 = barrettWithoutPrecomputation64_si128(
       _mm_clmulepi64_si128( tkey2, tkey2, 0x00));
    __m128i tkey4 = barrettWithoutPrecomputation64_si128(
       _mm_clmulepi64_si128( tkey3, tkey3, 0x00));
    // powers of the key are packed into two registers
    __m128i key = _mm_xor_si128(tkey1,_mm_slli_si128(tkey2,8));
    __m128i key2 = _mm_xor_si128(tkey3,_mm_slli_si128(tkey4,8));
    __m128i acc = _mm_set_epi64x(0,*string);
    __m128i mask = _mm_set_epi64x(0,-1);
    ++string;
    for(; string+3< endstring; string+=4 ) {
        __m128i temp = _mm_lddqu_si128((__m128i *) string);
        __m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2));
        const __m128i x1 =  _mm_and_si128 (temp,mask);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp, key, 0x01);
        const __m128i clprod2  = _mm_clmulepi64_si128( temp2, key, 0x10);
        const __m128i clprod3  = _mm_clmulepi64_si128( temp2, key2, 0x01);
        acc  = _mm_clmulepi64_si128( acc, key2, 0x10);
        acc = _mm_xor_si128 (acc,_mm_xor_si128 (_mm_xor_si128 (x1,clprod1),_mm_xor_si128 (clprod2,clprod3)));
    }
    for(; string+1< endstring; string+=2 ) {
        __m128i temp = _mm_lddqu_si128((__m128i *) string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp, key, 0x01);
        acc  = _mm_clmulepi64_si128( acc, key, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        acc = _mm_xor_si128 (acc,_mm_and_si128 (temp,mask));
        acc = barrettWithoutPrecomputation64_si128(acc);
    }
    if(string < endstring) {
        const __m128i temp = _mm_set_epi64x(0,*string);
        const __m128i multi = _mm_clmulepi64_si128( acc, key, 0x00);
        acc = barrettWithoutPrecomputation64_si128(multi);
        acc = _mm_xor_si128 (acc,temp);
    }
    __m128i multi = _mm_clmulepi64_si128( acc, key, 0x00);
    return barrettWithoutPrecomputation64(multi);
}

#endif

#endif /* CLMULHASHFUNCTIONS64BITS_H_ */
