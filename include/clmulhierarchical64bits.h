/*
 * clmulhierarchical64bits.h
 *
 *  Created on: Jul 29, 2014
 *      Author: lemire
 */

#ifndef CLMULHIERARCHICAL64BITS_H_
#define CLMULHIERARCHICAL64BITS_H_

#include "clmul.h"

///////////
// The main contribution of this header file is CLHASH
/////////

// multiply the length and the some key, no modulo
static __m128i lazyLengthHash(uint64_t keylength, uint64_t length) {
    const __m128i lengthvector = _mm_set_epi64x(keylength,length);
    const __m128i clprod1  = _mm_clmulepi64_si128( lengthvector, lengthvector, 0x10);
    return clprod1;
}
// hashing the bits in value using the keys key1 and key2 (only the first 64 bits of key2 are used).
// This is basically (a xor k1) * (b xor k2) mod p with length component
static uint64_t simple128to64hashwithlength( const __m128i value, const __m128i key, uint64_t keylength, uint64_t length) {
    const __m128i add =  _mm_xor_si128 (value,key);
    const __m128i clprod1  = _mm_clmulepi64_si128( add, add, 0x10);
    const __m128i total = _mm_xor_si128 (clprod1,lazyLengthHash(keylength, length));
    return precompReduction64(total);
}


enum {CLHASH_DEBUG=0};

// For use with CLHASH
// we expect length to have value 128 or, at least, to be divisible by 4.
static __m128i __clmulhalfscalarproductwithoutreduction(const __m128i * randomsource, const uint64_t * string,
        const size_t length) {
    assert(((uintptr_t) randomsource & 15) == 0);// we expect cache line alignment for the keys
    // we expect length = 128, so we need  16 cache lines of keys and 16 cache lines of strings.
    if(CLHASH_DEBUG) assert((length & 3) == 0); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length;
    __m128i acc = _mm_setzero_si128();
    // we expect length = 128
    for (; string + 3 < endstring; randomsource += 2, string += 4) {
        const __m128i temp1 = _mm_load_si128( randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 = _mm_xor_si128(temp1, temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
        acc = _mm_xor_si128(clprod1, acc);
        const __m128i temp12 = _mm_load_si128(randomsource + 1);
        const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string + 2));
        const __m128i add12 = _mm_xor_si128(temp12, temp22);
        const __m128i clprod12 = _mm_clmulepi64_si128(add12, add12, 0x10);
        acc = _mm_xor_si128(clprod12, acc);
    }
    if(CLHASH_DEBUG) assert(string == endstring);
    return acc;
}



// the value length does not have to be divisible by 4
static __m128i __clmulhalfscalarproductwithtailwithoutreduction(const __m128i * randomsource,
        const uint64_t * string, const size_t length) {
    assert(((uintptr_t) randomsource & 15) == 0);// we expect cache line alignment for the keys
    const uint64_t * const endstring = string + length;
    __m128i acc = _mm_setzero_si128();
    for (; string + 3 < endstring; randomsource += 2, string += 4) {
        const __m128i temp1 = _mm_load_si128(randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 = _mm_xor_si128(temp1, temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
        acc = _mm_xor_si128(clprod1, acc);
        const __m128i temp12 = _mm_load_si128(randomsource+1);
        const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string + 2));
        const __m128i add12 = _mm_xor_si128(temp12, temp22);
        const __m128i clprod12 = _mm_clmulepi64_si128(add12, add12, 0x10);
        acc = _mm_xor_si128(clprod12, acc);
    }
    if (string + 1 < endstring) {
        if(CLHASH_DEBUG) assert(length != 128);
        const __m128i temp1 = _mm_load_si128(randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 = _mm_xor_si128(temp1, temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
        acc = _mm_xor_si128(clprod1, acc);
        randomsource += 1;
        string += 2;
    }
    if (string < endstring) {
        if(CLHASH_DEBUG) assert(length != 128);
        const __m128i temp1 = _mm_load_si128(randomsource);
        const __m128i temp2 = _mm_loadl_epi64((__m128i const*)string);
        const __m128i add1 = _mm_xor_si128(temp1, temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
        acc = _mm_xor_si128(clprod1, acc);
        if(CLHASH_DEBUG) ++string;
    }
    if(CLHASH_DEBUG) assert(string == endstring);
    return acc;
}
// the value length does not have to be divisible by 4
static __m128i __clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(const __m128i * randomsource,
        const uint64_t * string, const size_t length, const uint64_t extraword) {
    assert(((uintptr_t) randomsource & 15) == 0);// we expect cache line alignment for the keys
    const uint64_t * const endstring = string + length;
    __m128i acc = _mm_setzero_si128();
    for (; string + 3 < endstring; randomsource += 2, string += 4) {
        const __m128i temp1 = _mm_load_si128(randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 = _mm_xor_si128(temp1, temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
        acc = _mm_xor_si128(clprod1, acc);
        const __m128i temp12 = _mm_load_si128(randomsource+1);
        const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string + 2));
        const __m128i add12 = _mm_xor_si128(temp12, temp22);
        const __m128i clprod12 = _mm_clmulepi64_si128(add12, add12, 0x10);
        acc = _mm_xor_si128(clprod12, acc);
    }
    if (string + 1 < endstring) {
        const __m128i temp1 = _mm_load_si128(randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 = _mm_xor_si128(temp1, temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
        acc = _mm_xor_si128(clprod1, acc);
        randomsource += 1;
        string += 2;
    }
    // we have to append an extra 1
    if (string < endstring) {
        const __m128i temp1 = _mm_load_si128(randomsource);
        const __m128i temp2 = _mm_set_epi64x(extraword,*string);
        const __m128i add1 = _mm_xor_si128(temp1, temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
        acc = _mm_xor_si128(clprod1, acc);
    } else {
        const __m128i temp1 = _mm_load_si128(randomsource);
        const __m128i temp2 = _mm_loadl_epi64((__m128i const*)&extraword);
        const __m128i add1 = _mm_xor_si128(temp1, temp2);
        const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x01);
        acc = _mm_xor_si128(clprod1, acc);
    }
    return acc;
}


static __m128i __clmulhalfscalarproductOnlyExtraWord(const __m128i * randomsource,
        const uint64_t extraword) {
    const __m128i temp1 = _mm_load_si128(randomsource);
    const __m128i temp2 = _mm_loadl_epi64((__m128i const*)&extraword);
    const __m128i add1 = _mm_xor_si128(temp1, temp2);
    const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x01);
    return clprod1;
}




#ifdef BITMIX
////////
// an invertible function used to mix the bits
// borrowed directly from murmurhash
////////
inline uint64_t fmix64 ( uint64_t k ) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}
#endif

enum {RANDOM_64BITWORDS_NEEDED_FOR_CLHASH=133,RANDOM_BYTES_NEEDED_FOR_CLHASH=133*8};

//////////////////////
// just two levels like VHASH
// at low level, we use a half-multiplication multilinear that we aggregate using
// a CLMUL polynomial hash
// this (RANDOM_BYTES_NEEDED_FOR_CLHASH random bytes or about 1KB)
//
// rs : the random data source (should contain at least RANDOM_BYTES_NEEDED_FOR_CLHASH random bytes)
// string : the input data source
// length : number of 64-bit words in the string
//////////////////////
uint64_t CLHASH(const void* rs, const uint64_t * string,
                const size_t length) {
    assert(sizeof(size_t)<=sizeof(uint64_t));// otherwise, we need to worry
    assert(((uintptr_t) rs & 15) == 0);// we expect cache line alignment for the keys
    const unsigned int m = 128;// we process the data in chunks of 16 cache lines
    if(CLHASH_DEBUG) assert((m  & 3) == 0); //m should be divisible by 4
    const int m128neededperblock = m / 2;// that is how many 128-bit words of random bits we use per block
    const __m128i * rs64 = (__m128i *) rs;
    __m128i polyvalue =  _mm_load_si128(rs64 + m128neededperblock); // to preserve alignment on cache lines for main loop, we pick random bits at the end
    polyvalue = _mm_and_si128(polyvalue,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x3fffffff));// setting two highest bits to zero
    // we should check that polyvalue is non-zero, though this is best done outside the function and highly unlikely
    if (m < length) { // long strings
        __m128i  acc =  __clmulhalfscalarproductwithoutreduction(rs64, string,m);
        size_t t = m;
        for (; t +  m <= length; t +=  m) {
            // we compute something like
            // acc+= polyvalue * acc + h1
            acc =  mul128by128to128_lazymod127(polyvalue,acc);
            const __m128i h1 =  __clmulhalfscalarproductwithoutreduction(rs64, string+t,m);
            acc = _mm_xor_si128(acc,h1);
        }
        const int remain = length - t;
        // we compute something like
        // acc+= polyvalue * acc + h1
        acc = mul128by128to128_lazymod127(polyvalue, acc);
        const __m128i h1 =
            __clmulhalfscalarproductwithtailwithoutreduction(
                rs64, string + t, remain);
        acc = _mm_xor_si128(acc, h1);
        const __m128i finalkey = _mm_load_si128(rs64 + m128neededperblock + 1);
        const uint64_t keylength = *(const uint64_t *)(rs64 + m128neededperblock + 2);
        return simple128to64hashwithlength(acc,finalkey,keylength, (uint64_t)(length * sizeof(uint64_t)));
    } else { // short strings
        __m128i  acc = __clmulhalfscalarproductwithtailwithoutreduction(rs64, string, length);
        const uint64_t keylength = *(const uint64_t *)(rs64 + m128neededperblock + 2);
        acc = _mm_xor_si128(acc,lazyLengthHash(keylength, (uint64_t)(length * sizeof(uint64_t))));
#ifdef BITMIX
        return fmix64(precompReduction64(acc)) ;
#else
        return precompReduction64(acc) ;
#endif
    }
}
// there always remain an incomplete word that has 1,2, 3, 4, 5, 6, 7 used bytes.
// we append 0s to it
static inline uint64_t createLastWord(const size_t lengthbyte, const uint64_t * lastw) {
    const int significantbytes = lengthbyte % sizeof(uint64_t);
    uint64_t lastword = 0;
    memcpy(&lastword,lastw,significantbytes); // could possibly be faster?
    return lastword;
}

//////////////////////
// like CLHASH, but can hash byte strings
//
// rs : the random data source (should contain at least RANDOM_BYTES_NEEDED_FOR_CLHASH random bytes)
// stringbyte : the input data source
// length : number of bytes in the string
//////////////////////
uint64_t CLHASHbyte(const void* rs, const char * stringbyte,
                    const size_t lengthbyte) {
    assert(sizeof(size_t)<=sizeof(uint64_t));// otherwise, we need to worry
    assert(((uintptr_t) rs & 15) == 0);// we expect cache line alignment for the keys
    const unsigned int  m = 128;// we process the data in chunks of 16 cache lines
    if(CLHASH_DEBUG) assert((m  & 3) == 0); //m should be divisible by 4
    const int m128neededperblock = m / 2;// that is how many 128-bit words of random bits we use per block
    const __m128i * rs64 = (__m128i *) rs;
    __m128i polyvalue =  _mm_load_si128(rs64 + m128neededperblock); // to preserve alignment on cache lines for main loop, we pick random bits at the end
    polyvalue = _mm_and_si128(polyvalue,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x3fffffff));// setting two highest bits to zero
    // we should check that polyvalue is non-zero, though this is best done outside the function and highly unlikely
    const size_t length = lengthbyte / sizeof(uint64_t); // # of complete words
    const size_t lengthinc = (lengthbyte + sizeof(uint64_t) - 1) / sizeof(uint64_t); // # of words, including partial ones

    const uint64_t * string = (const uint64_t *)  stringbyte;
    if (m < lengthinc) { // long strings // modified from length to lengthinc to address issue #3 raised by Eik List
        __m128i  acc =  __clmulhalfscalarproductwithoutreduction(rs64, string,m);
        size_t t = m;
        for (; t +  m <= length; t +=  m) {
            // we compute something like
            // acc+= polyvalue * acc + h1
            acc =  mul128by128to128_lazymod127(polyvalue,acc);
            const __m128i h1 =  __clmulhalfscalarproductwithoutreduction(rs64, string+t,m);
            acc = _mm_xor_si128(acc,h1);
        }
        const int remain = length - t;  // number of completely filled words

        if (remain != 0) {
            // we compute something like
            // acc+= polyvalue * acc + h1
            acc = mul128by128to128_lazymod127(polyvalue, acc);
            if (lengthbyte % sizeof(uint64_t) == 0) {
                const __m128i h1 =
                    __clmulhalfscalarproductwithtailwithoutreduction(rs64,
                            string + t, remain);
                acc = _mm_xor_si128(acc, h1);
            } else {
                const uint64_t lastword = createLastWord(lengthbyte,
                                          (string + length));
                const __m128i h1 =
                    __clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(
                        rs64, string + t, remain, lastword);
                acc = _mm_xor_si128(acc, h1);
            }
        } else if (lengthbyte % sizeof(uint64_t) != 0) {// added to address issue #2 raised by Eik List
            // there are no completely filled words left, but there is one partial word.
            acc = mul128by128to128_lazymod127(polyvalue, acc);
            const uint64_t lastword = createLastWord(lengthbyte, (string + length));
            const __m128i h1 = __clmulhalfscalarproductOnlyExtraWord( rs64, lastword);
            acc = _mm_xor_si128(acc, h1);
        }

        const __m128i finalkey = _mm_load_si128(rs64 + m128neededperblock + 1);
        const uint64_t keylength = *(const uint64_t *)(rs64 + m128neededperblock + 2);
        return simple128to64hashwithlength(acc,finalkey,keylength, (uint64_t)lengthbyte);
    } else { // short strings
        if(lengthbyte % sizeof(uint64_t) == 0) {
            __m128i  acc = __clmulhalfscalarproductwithtailwithoutreduction(rs64, string, length);
            const uint64_t keylength = *(const uint64_t *)(rs64 + m128neededperblock + 2);
            acc = _mm_xor_si128(acc,lazyLengthHash(keylength, (uint64_t)lengthbyte));
#ifdef BITMIX
            return fmix64(precompReduction64(acc)) ;
#else
            return precompReduction64(acc) ;
#endif
        }
        const uint64_t lastword = createLastWord(lengthbyte, (string + length));
        __m128i acc = __clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(
                          rs64, string, length, lastword);
        const uint64_t keylength =  *(const uint64_t *)(rs64 + m128neededperblock + 2);
        acc = _mm_xor_si128(acc,lazyLengthHash(keylength, (uint64_t)lengthbyte));
#ifdef BITMIX
        return fmix64(precompReduction64(acc)) ;
#else
        return precompReduction64(acc) ;
#endif
    }
}


#endif /* CLMULHIERARCHICAL64BITS_H_ */
