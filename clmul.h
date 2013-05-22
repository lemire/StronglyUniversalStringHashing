#ifndef _CLMUL_H_
#define _CLMUL_H_


#ifdef __PCLMUL__

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <wmmintrin.h>


/////////////////////////////////////////////////////////////////
// working from 
// "Modular Reduction in GF(2n) without Pre-computational Phase"
// by M. Knezevic, K. Sakiyama, J. Fan, and I. Verbauwhede (2008)
// algo. 4
//
// A is input, M is irred. (n=32)
//
//(((A div x^n) * M ) div x^n) * M) mod x^n
//+(A mod x^n)
//////////////////////////////////////////////////
uint32_t barrettWithoutPrecomputation( __m128i A) {
	///http://www.jjj.de/mathdata/minweight-primpoly.txt
	const uint64_t irredpoly = 1UL+(1UL<<2)+(1UL<<6)+(1UL<<7)+(1UL<<32);
	// it is important, for the algo. we have chosen that 7 is smaller 
	// equal than 16=32/2
	const int n = 32;// degree of the polynomial
	const __m128i C = _mm_set_epi64x(0,irredpoly);// C is the irreducible poly.
	/////////////////
	/// This algo. requires two multiplications (_mm_clmulepi64_si128)
	/// They are probably the bottleneck.
	/// Note: Barrett's original algorithm also required two multiplications.
	////////////////
	__m128i Q1 = _mm_srli_epi64 (A, n);
	__m128i Q2 = _mm_clmulepi64_si128( Q1, C, 0x00);// A div x^n
	__m128i Q3 = _mm_srli_epi64 (Q2, n);
	// commenting out the long way derived from the paper (following two lines are enough)
	//__m128i R1 = _mm_and_si128 (maskm128,A);
	//__m128i R2 = _mm_and_si128 (maskm128,_mm_clmulepi64_si128( Q3, C, 0x00));
	//__m128i final  = _mm_xor_si128 (R1, R2);
	__m128i final  = _mm_xor_si128 (A, _mm_clmulepi64_si128( Q3, C, 0x00));		
	return _mm_cvtsi128_si64(final); 
}


uint32_t hashGaloisFieldMultilinear(const uint64_t *  randomsource, const uint32_t *  string, const size_t length) {
	const uint32_t * const endstring = string + length;
	const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    __m128i acc = _mm_set_epi64x(0,*(randomsource32++));
    for(; string!= endstring; ++randomsource32,++string ) {
    	__m128i temp = _mm_set_epi64x(*randomsource32,*string);
    	__m128i clprod  = _mm_clmulepi64_si128( temp, temp, 0x10);
        acc = _mm_xor_si128 (clprod,acc);   
    }
    return barrettWithoutPrecomputation(acc);
}

uint32_t hashGaloisFieldMultilinearHalfMultiplications(const uint64_t*  randomsource, const uint32_t *  string, const size_t length) {
	assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
	const uint32_t * const endstring = string + length;
	const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    __m128i acc = _mm_set_epi64x(0,*(randomsource32++));
    for(; string!= endstring; randomsource32+=2,string+=2 ) {
    	__m128i temp1 = _mm_set_epi64x(*randomsource32,*(randomsource32+1));
    	__m128i temp2 = _mm_set_epi64x(*string,*(string+1));
    	__m128i twosums = _mm_xor_si128(temp1,temp2); 
    	__m128i clprod  = _mm_clmulepi64_si128( twosums, twosums, 0x10);
        acc = _mm_xor_si128 (clprod,acc);   
    }
    return barrettWithoutPrecomputation(acc);
}

uint32_t hashGaloisFieldfast(const uint64_t*  randomsource, const uint32_t *  string, const size_t length) {
	assert(length / 4 * 4 == length); // if not, we need special handling (omitted)
	const uint32_t * const endstring = string + length;
	const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    __m128i acc = _mm_set_epi64x(0,*(randomsource32++));
    const __m128i zero =  _mm_setzero_si128 ();
    for(; string!= endstring; randomsource32+=4,string+=4 ) {
    	const __m128i temp1 = _mm_load_si128((__m128i * )randomsource32);
    	const __m128i temp2 = _mm_load_si128((__m128i *) string);
    	const __m128i twosums = _mm_xor_si128(temp1,temp2); 
    	const __m128i part1 = _mm_unpacklo_epi32(twosums,zero);
	const __m128i clprod1  = _mm_clmulepi64_si128( part1, part1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);   
    	const __m128i part2 = _mm_unpackhi_epi32(twosums,zero);
	const __m128i clprod2  = _mm_clmulepi64_si128( part2, part2, 0x10);
        acc = _mm_xor_si128 (clprod2,acc);   
     }
    return barrettWithoutPrecomputation(acc);
}

#endif

#endif
