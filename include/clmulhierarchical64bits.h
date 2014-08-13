/*
 * clmulhierarchical64bits.h
 *
 *  Created on: Jul 29, 2014
 *      Author: lemire
 */

#ifndef CLMULHIERARCHICAL64BITS_H_
#define CLMULHIERARCHICAL64BITS_H_



// For use with hashCLMUL2Level
__m128i __clmulhalfscalarproductwithoutreduction(const void* rs, const uint64_t * string,
		const size_t length) {
	assert(length / 4 * 4 == length); // if not, we need special handling (omitted)
	const uint64_t * const endstring = string + length;
	const uint64_t * randomsource = (const uint64_t *) rs;
	__m128i acc = _mm_set_epi64x(0, *(randomsource));
	randomsource += 1;
	for (; string + 3 < endstring; randomsource += 4, string += 4) {
		const __m128i temp1 = _mm_lddqu_si128((__m128i *) randomsource);
		const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
		const __m128i add1 = _mm_xor_si128(temp1, temp2);
		const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
		acc = _mm_xor_si128(clprod1, acc);
		const __m128i temp12 = _mm_lddqu_si128((__m128i *) (randomsource + 2));
		const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string + 2));
		const __m128i add12 = _mm_xor_si128(temp12, temp22);
		const __m128i clprod12 = _mm_clmulepi64_si128(add12, add12, 0x10);
		acc = _mm_xor_si128(clprod12, acc);
	}
	return acc;
}


// For use with hashCLMUL2Level
__m128i __clmulhalfscalarproductwithtailwithoutreduction(const void* rs,
		const uint64_t * string, const size_t length) {
	const uint64_t * const endstring = string + length;
	const uint64_t * randomsource = (const uint64_t *) rs;
	__m128i acc = _mm_set_epi64x(0, *(randomsource));
	randomsource += 1;
	for (; string + 3 < endstring; randomsource += 4, string += 4) {
		const __m128i temp1 = _mm_lddqu_si128((__m128i *) randomsource);
		const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
		const __m128i add1 = _mm_xor_si128(temp1, temp2);
		const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
		acc = _mm_xor_si128(clprod1, acc);
		const __m128i temp12 = _mm_lddqu_si128((__m128i *) (randomsource + 2));
		const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string + 2));
		const __m128i add12 = _mm_xor_si128(temp12, temp22);
		const __m128i clprod12 = _mm_clmulepi64_si128(add12, add12, 0x10);
		acc = _mm_xor_si128(clprod12, acc);
	}
	if (string + 1 < endstring) {
		const __m128i temp1 = _mm_lddqu_si128((__m128i *) randomsource);
		const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
		const __m128i add1 = _mm_xor_si128(temp1, temp2);
		const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
		acc = _mm_xor_si128(clprod1, acc);
		randomsource += 2;
		string += 2;
	}
	if (string < endstring) {
		const __m128i temp1 = _mm_set_epi64x(0, *randomsource);
		const __m128i temp2 = _mm_set_epi64x(0, *string);
		const __m128i clprod1 = _mm_clmulepi64_si128(temp1, temp2, 0x00);
		acc = _mm_xor_si128(clprod1, acc);
	}
	return acc;
}



// just two levels like VHASH
// at low level, we use multilinear that we aggregate using
// a CLMUL polynomial hash
// this uses 128 + 1 keys.(129*8 random bytes or about 1KB)
// *WARNING*  this is highly experimental TODO: validate the result
uint64_t hashCLMUL2Level(const void* rs, const uint64_t * string,
		const size_t length) {
	if (length == 0)
		return 0; // hmmmm...
	const int m = 128;
	__m128i polyvalue =  _mm_lddqu_si128((__m128i *) rs);
	polyvalue = _mm_and_si128(polyvalue,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x3fffffff));// setting two highest bits to zero
	// we should check that polyvalue is non-zero
	const __m128i * rs64 = (__m128i *) rs + 1;

	// TODO: To handle variable length inputs, we should XOR the result with
	// the hash value of the length

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
		//acc = mod127fromlazy(acc);
		return _mm_cvtsi128_si64(acc);// TODO: This is not correct, we need to hash it down (like VHASH)
	} else { // short strings
		__m128i  acc = __clmulhalfscalarproductwithtailwithoutreduction(rs64, string, length);
		//acc = mod127fromlazy(acc);
		return _mm_cvtsi128_si64(acc);// TODO: This is not correct, we need to hash it down (like VHASH)
	}
}

#endif /* CLMULHIERARCHICAL64BITS_H_ */
