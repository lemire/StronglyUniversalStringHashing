/*
 * clmulhierarchical64bits.h
 *
 *  Created on: Jul 29, 2014
 *      Author: lemire
 */

#ifndef CLMULHIERARCHICAL64BITS_H_
#define CLMULHIERARCHICAL64BITS_H_

#include "clmul.h"
//#include "iacaMarks.h"

///////////
// The main contribution of this header file is CLHASH
/////////

// hashing the bits in value using the keys key1 and key2 (only the first 64 bits of key2 are used).
// This is basically (a xor k1) * (b xor k2) mod p
//
static uint64_t simple128to64hash( __m128i value, __m128i key) {
    const __m128i add =  _mm_xor_si128 (value,key);
    const __m128i clprod1  = _mm_clmulepi64_si128( add, add, 0x10);
	return precompReduction64(clprod1);
}


// For use with CLHASH
// we expect length to have value 128 or, at least, to be divisible by 4.
static __m128i __clmulhalfscalarproductwithoutreduction(const __m128i * randomsource, const uint64_t * string,
		const size_t length) {
	assert(((uintptr_t) randomsource & 15) == 0);// we expect cache line alignment for the keys
	// we expect length = 128, so we need  16 cache lines of keys and 16 cache lines of strings.
	assert((length & 3) == 0); // if not, we need special handling (omitted)
	const uint64_t * const endstring = string + length;
	__m128i acc = _mm_setzero_si128();
	// we expect length = 128
	for (; string + 3 < endstring; randomsource += 2, string += 4) {
          //IACA_START
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
        //IACA_END
	assert(string == endstring);
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
		assert(length != 128);
		const __m128i temp1 = _mm_load_si128(randomsource);
		const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
		const __m128i add1 = _mm_xor_si128(temp1, temp2);
		const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
		acc = _mm_xor_si128(clprod1, acc);
		randomsource += 1;
		string += 2;
	}
	if (string < endstring) {
		assert(length != 128);
		const __m128i temp1 = _mm_load_si128(randomsource);
		const __m128i temp2 = _mm_loadl_epi64((__m128i const*)string);
		const __m128i add1 = _mm_xor_si128(temp1, temp2);
		const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
		acc = _mm_xor_si128(clprod1, acc);
		++string;
	}
	assert(string == endstring);
	return acc;
}
// the value length does not have to be divisible by 4
static __m128i __clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(const __m128i * randomsource,
		const uint64_t * string, const size_t length,uint64_t extraword) {
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


//////////////////////
// just two levels like VHASH
// at low level, we use a half-multiplication multilinear that we aggregate using
// a CLMUL polynomial hash
// this uses 128 + 2 keys.(RANDOM_BYTES_NEEDED_FOR_CLHASH random bytes or about 1KB)
//
// rs : the random data source (should contain at least 130*8 random bytes)
// string : the input data source
// length : number of 64-bit words in the string
//////////////////////
uint64_t CLHASH(const void* rs, const uint64_t * string,
		const size_t length) {
	assert(sizeof(size_t)<=sizeof(uint64_t));// otherwise, we need to worry
	assert(((uintptr_t) rs & 15) == 0);// we expect cache line alignment for the keys
	const unsigned int m = 128;// we process the data in chunks of 16 cache lines
	assert((m  & 3) == 0); //m should be divisible by 4
	const int m128neededperblock = m / 2;// that is how many 128-bit words of random bits we use per block
	const __m128i * rs64 = (__m128i *) rs;
	__m128i polyvalue =  _mm_load_si128(rs64 + m128neededperblock); // to preserve alignment on cache lines for main loop, we pick random bits at the end
	polyvalue = _mm_and_si128(polyvalue,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x3fffffff));// setting two highest bits to zero
	// we should check that polyvalue is non-zero, though this is best done outside the function and highly unlikely
	if (m <= length) { // long strings
		__m128i  acc =  __clmulhalfscalarproductwithoutreduction(rs64, string,m);
		size_t t = m;
		for (; t +  m <= length; t +=  m) {
			// we compute something like
			// acc+= polyvalue * acc + h1
			acc =  mul128by128to128_lazymod127(polyvalue,acc);
			__m128i h1 =  __clmulhalfscalarproductwithoutreduction(rs64, string+t,m);
			acc = _mm_xor_si128(acc,h1);
		}
		int remain = length - t;
		{
			// we compute something like
			// acc+= polyvalue * acc + h1
			acc =  mul128by128to128_lazymod127(polyvalue,acc);
			uint64_t lastword = 1;
			__m128i h1 =  __clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(rs64, string+t,remain,lastword);
			acc = _mm_xor_si128(acc,h1);
		}
		__m128i finalkey = _mm_load_si128(rs64 + m128neededperblock + 1);
		return  simple128to64hash(acc,finalkey );
	} else { // short strings
		uint64_t lastword = 1;
		__m128i  acc = __clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(rs64, string, length, lastword);
		return fmix64(precompReduction64(acc)) ;
		//fmix64(
	}
}


// there always remain an incomplete word that has 0,1,2, 3, 4, 5, 6, 7 used bytes.
// we append 1 to it
uint64_t createLastWord(size_t lengthbyte, uint64_t lastw) {
	int significantbytes = lengthbyte % sizeof(uint64_t);
	if(significantbytes==0) return 1;
	uint64_t mask = (~((uint64_t)0)) >> ((sizeof(uint64_t)- significantbytes)*8);
	uint64_t lastword = lastw  & mask;// could be cleverer
	lastword |= ((uint64_t)1) <<( significantbytes  * 8);
	return lastword;
}

// there always remain an incomplete word that has 0,1,2, 3, 4, 5, 6, 7 used bytes.
// we append 1 to it
uint64_t createUnpaddedLastWord(size_t lengthbyte, uint64_t lastw) {
	int significantbytes = lengthbyte % sizeof(uint64_t);
	if(significantbytes==0) return 0;
	uint64_t mask = (~((uint64_t)0)) >> ((sizeof(uint64_t)- significantbytes)*8);
	uint64_t lastword = lastw  & mask;// could be cleverer
	return lastword;
}

enum{RANDOM_BYTES_NEEDED_FOR_CLHASH=132*8};

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
	assert((m  & 3) == 0); //m should be divisible by 4
	const int m128neededperblock = m / 2;// that is how many 128-bit words of random bits we use per block
	const __m128i * rs64 = (__m128i *) rs;
	__m128i polyvalue =  _mm_load_si128(rs64 + m128neededperblock); // to preserve alignment on cache lines for main loop, we pick random bits at the end
	polyvalue = _mm_and_si128(polyvalue,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x3fffffff));// setting two highest bits to zero
	// we should check that polyvalue is non-zero, though this is best done outside the function and highly unlikely
	size_t length = lengthbyte / sizeof(uint64_t);
	const uint64_t * string = (const uint64_t *)  stringbyte;
	if (m <= length) { // long strings
		__m128i  acc =  __clmulhalfscalarproductwithoutreduction(rs64, string,m);
		size_t t = m;
		for (; t +  m <= length; t +=  m) {
			// we compute something like
			// acc+= polyvalue * acc + h1
			acc =  mul128by128to128_lazymod127(polyvalue,acc);
			__m128i h1 =  __clmulhalfscalarproductwithoutreduction(rs64, string+t,m);
			acc = _mm_xor_si128(acc,h1);
		}
		int remain = length - t;
		{
			// we compute something like
			// acc+= polyvalue * acc + h1
			acc =  mul128by128to128_lazymod127(polyvalue,acc);
			uint64_t lastword = createLastWord(lengthbyte, * (string + length));
			__m128i h1 =  __clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(rs64, string+t,remain,lastword);
			acc = _mm_xor_si128(acc,h1);
		}
		__m128i finalkey = _mm_load_si128(rs64 + m128neededperblock + 1);
		return  simple128to64hash(acc,finalkey );
	} else { // short strings
		uint64_t lastword = createLastWord(lengthbyte, * (string + length));
		__m128i  acc = __clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(rs64, string, length, lastword);
		return  fmix64(precompReduction64(acc)) ;
	}
}

//////////////////////
// like CLHASHbyte but does not pad with an extra 1
//
// rs : the random data source (should contain at least RANDOM_BYTES_NEEDED_FOR_CLHASH random bytes)
// stringbyte : the input data source
// length : number of bytes in the string
//////////////////////
uint64_t CLHASHbyteFixed(const void* rs, const char * stringbyte,
		const size_t lengthbyte) {
	assert(sizeof(size_t) <= sizeof(uint64_t)); // otherwise, we need to worry
	assert(((uintptr_t) rs & 15) == 0); // we expect cache line alignment for the keys
	const unsigned int m = 128; // we process the data in chunks of 16 cache lines
	assert((m & 3) == 0); //m should be divisible by 4
	const int m128neededperblock = m / 2; // that is how many 128-bit words of random bits we use per block
	const __m128i * rs64 = (__m128i *) rs;
	__m128i polyvalue = _mm_load_si128(rs64 + m128neededperblock); // to preserve alignment on cache lines for main loop, we pick random bits at the end
	polyvalue = _mm_and_si128(polyvalue,
			_mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x3fffffff)); // setting two highest bits to zero
	// we should check that polyvalue is non-zero, though this is best done outside the function and highly unlikely
	size_t length = lengthbyte / sizeof(uint64_t);
	const uint64_t * string = (const uint64_t *) stringbyte;
	if (m < length) { // long strings
		__m128i acc = __clmulhalfscalarproductwithoutreduction(rs64, string, m);
		size_t t = m;
		for (; t + m <= length; t += m) {
			// we compute something like
			// acc+= polyvalue * acc + h1
			acc = mul128by128to128_lazymod127(polyvalue, acc);
			__m128i h1 = __clmulhalfscalarproductwithoutreduction(rs64,
					string + t, m);
			acc = _mm_xor_si128(acc, h1);
		}
		int remain = length - t;
		if(remain != 0){
			// we compute something like
			// acc+= polyvalue * acc + h1
			acc = mul128by128to128_lazymod127(polyvalue, acc);
			if (lengthbyte % sizeof(uint64_t) == 0) {
				__m128i h1 = __clmulhalfscalarproductwithtailwithoutreduction(
						rs64, string + t, remain);
				acc = _mm_xor_si128(acc, h1);

			} else {

				uint64_t lastword = createUnpaddedLastWord(lengthbyte,
						*(string + length));
				__m128i h1 =
						__clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(
								rs64, string + t, remain, lastword);
				acc = _mm_xor_si128(acc, h1);
			}
		}
		__m128i finalkey = _mm_load_si128(rs64 + m128neededperblock + 1);
		return simple128to64hash(acc, finalkey);
	} else { // short strings
		if(lengthbyte % sizeof(uint64_t) == 0) {
			__m128i acc = __clmulhalfscalarproductwithtailwithoutreduction(
									rs64, string, length);
			return  fmix64(precompReduction64(acc)) ;
		}
		uint64_t lastword = createUnpaddedLastWord(lengthbyte,
				*(string + length));
		__m128i acc = __clmulhalfscalarproductwithtailwithoutreductionWithExtraWord(
						rs64, string, length, lastword);
		return  fmix64(precompReduction64(acc)) ;
	}
}

/////////
// what follows are convenience functions
// call init_clhash once with a 32-bit key
// then call clhash to hash strings.
/////////
#include "mersenne.h"

static uint64_t randomkey[RANDOM_BYTES_NEEDED_FOR_CLHASH];

void init_clhash( uint32_t seed) {
  ZRandom zr;
  initZRandom(&zr,seed);
  for (int i=0; i < RANDOM_BYTES_NEEDED_FOR_CLHASH; ++i)
	  randomkey[i] = getValue(&zr) | ( ((uint64_t) getValue(&zr)) << 32);
}

uint64_t clhash( const void *key, int len) {
	return CLHASHbyte(randomkey,(const char *)key,len);
}

#endif /* CLMULHIERARCHICAL64BITS_H_ */
