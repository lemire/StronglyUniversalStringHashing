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
// The main contribution of this header file is hashCLMUL2Level
/////////

// hashing the bits in value using the keys key1 and key2 (only the first 64 bits of key2 are used).
// This is basically (a xor k1) * (b xor k2) mod p
//
static uint64_t simple128to64hash( __m128i value, __m128i key) {
    const __m128i add =  _mm_xor_si128 (value,key);
    const __m128i clprod1  = _mm_clmulepi64_si128( add, add, 0x10);
	return precompReduction64(clprod1);
}


// For use with hashCLMUL2Level
// we expect length to have value 128 or, at least, to be divisible by 4.
static __m128i __clmulhalfscalarproductwithoutreduction(const __m128i * randomsource, const uint64_t * string,
		const size_t length) {
	assert(((uintptr_t) randomsource & 15) == 0);// we expect cache line alignment for the keys
	// we expect length = 128, so we need  16 cache lines of keys and 16 cache lines of strings.
	assert((length & 3) == 0); // if not, we need special handling (omitted)
	const uint64_t * const endstring = string + length;
	__m128i acc = _mm_setzero_si128();


	// we process we expect length = 128
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
	return acc;
}


// For use with hashCLMUL2Level
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
		const __m128i temp1 = _mm_load_si128(randomsource);
		const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
		const __m128i add1 = _mm_xor_si128(temp1, temp2);
		const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
		acc = _mm_xor_si128(clprod1, acc);
		randomsource += 1;
		string += 2;
	}
	if (string < endstring) {
		const __m128i temp1 = _mm_load_si128(randomsource);
		const __m128i temp2 = _mm_loadl_epi64((__m128i const*)string);
		const __m128i clprod1 = _mm_clmulepi64_si128(temp1, temp2, 0x00);
		acc = _mm_xor_si128(clprod1, acc);
	}
	return acc;
}



// just two levels like VHASH
// at low level, we use a half-multiplication multilinear that we aggregate using
// a CLMUL polynomial hash
// this uses 128 + 1 keys.(129*8 random bytes or about 1KB)
uint64_t hashCLMUL2Level(const void* rs, const uint64_t * string,
		const size_t length) {
	if (length == 0)
		return 0; // hmmmm...
	assert(sizeof(size_t)<=sizeof(uint64_t));// otherwise, we need to worry
	assert(((uintptr_t) rs & 15) == 0);// we expect cache line alignment for the keys
	const int m = 128;// we process the data in chunks of 16 cache lines
	assert((m  & 3) == 0); //m should be divisible by 4
	const int m128neededperblock = m / 2;// that is how many 128-bit words of random bits we use per block
	const __m128i * rs64 = (__m128i *) rs;
	__m128i polyvalue =  _mm_load_si128(rs64 + m128neededperblock); // to preserve alignment on cache lines for main loop, we pick random bits at the end
	polyvalue = _mm_and_si128(polyvalue,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x3fffffff));// setting two highest bits to zero
	// we should check that polyvalue is non-zero, though this is best done outside the function and highly unlikely
	if (m < length) { // long strings
		__m128i  acc =  __clmulhalfscalarproductwithoutreduction(rs64, string,m);
		size_t t = m;
		for (; t +  m <= length; t +=  m) {
			acc =  mul128by128to128_lazymod127(polyvalue,acc);
			__m128i h1 =  __clmulhalfscalarproductwithoutreduction(rs64, string+t,m);
			acc = _mm_and_si128(acc,h1);
		}
		int remain = length - t;
		if(remain > 0) {
			acc =  mul128by128to128_lazymod127(polyvalue,acc);
			__m128i h1 =  __clmulhalfscalarproductwithtailwithoutreduction(rs64, string+t,remain);
			acc = _mm_and_si128(acc,h1);
		}
		return length ^ simple128to64hash(acc, _mm_load_si128(rs64 + m128neededperblock + 1));
	} else { // short strings
		__m128i  acc = __clmulhalfscalarproductwithtailwithoutreduction(rs64, string, length);
		return length ^ simple128to64hash(acc, _mm_load_si128(rs64 + m128neededperblock + 1));
	}
}


#endif /* CLMULHIERARCHICAL64BITS_H_ */
