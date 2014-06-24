#ifndef _CLMUL_H_
#define _CLMUL_H_
#ifdef __AVX__ // intel does not define PCLMUL
#define __PCLMUL__ 1
#endif

#ifdef __PCLMUL__

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <x86intrin.h>

//#define IACA
#ifdef IACA
#include </opt/intel/iaca/include/iacaMarks.h>
#else
#define IACA_START
#define IACA_END
#endif





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
/// WARNING: HIGH 96 BITS CONTAIN GARBAGE, must call _mm_cvtsi128_si32 to get
/// meaningful bits.
//// It assumes that only the lesser 64 bits are used.
__m128i barrettWithoutPrecomputation32_si128( __m128i A) {
	///http://www.jjj.de/mathdata/minweight-primpoly.txt
    const uint64_t irredpoly = 1UL+(1UL<<2)+(1UL<<6)+(1UL<<7)+(1UL<<32);
    // it is important, for the algo. we have chosen that 7 is smaller
    // equal than 16=32/2
    const int n = 32;// degree of the polynomial
    const __m128i C = _mm_set_epi64x(0,irredpoly);// C is the irreducible poly.
    /////////////////
    /// This algo. requires two multiplications (_mm_clmulepi64_si128)
    /// They are probably the bottleneck.
    /// Note: Barrett's original algorithm also requires two multiplications.
    ////////////////
    const __m128i Q1 = _mm_srli_si128 (A, 4);
    const __m128i Q2 = _mm_clmulepi64_si128( Q1, C, 0x00);// A div x^n
    const __m128i Q3 = _mm_srli_si128 (Q2, 4);
    const __m128i Q4 = _mm_clmulepi64_si128( Q3, C, 0x00);
    const __m128i final  = _mm_xor_si128 (A, Q4);
    return final; /// WARNING: HIGH 96 BITS CONTAIN GARBAGE
}

uint32_t barrettWithoutPrecomputation32( __m128i A) {
	return _mm_cvtsi128_si32(barrettWithoutPrecomputation32_si128(A));
}








/// WARNING: HIGH 64 BITS CONTAIN GARBAGE, must call _mm_cvtsi128_si64 to get
/// meaningful bits.
__m128i barrettWithoutPrecomputation64_si128( __m128i A) {
     ///http://www.jjj.de/mathdata/minweight-primpoly.txt
    // it is important, for the algo. we have chosen that 4 is smaller
    // equal than 32=64/2
    const int n = 64;// degree of the polynomial
    const __m128i C = _mm_set_epi64x(1U,(1U<<4)+(1U<<3)+(1U<<1)+(1U<<0));// C is the irreducible poly. (64,4,3,1,0)
    /////////////////
    /// This algo. requires two multiplications (_mm_clmulepi64_si128)
    /// They are probably the bottleneck.
    /// Note: Barrett's original algorithm also required two multiplications.
    ////////////////
    assert(n/8==8);
    __m128i Q2 = _mm_clmulepi64_si128( A, C, 0x01);
    Q2 = _mm_xor_si128(Q2,A);
    const __m128i Q4 = _mm_clmulepi64_si128( Q2, C, 0x01);
    const __m128i final  = _mm_xor_si128 (A, Q4);
    return final;/// WARNING: HIGH 64 BITS CONTAIN GARBAGE
 }
uint64_t barrettWithoutPrecomputation64( __m128i A) {
    const __m128i final  = barrettWithoutPrecomputation64_si128(A);
    return _mm_cvtsi128_si64(final);
 }




__m128i precompReduction64_si128( __m128i A) {
    const __m128i C = _mm_set_epi64x(1U,(1U<<4)+(1U<<3)+(1U<<1)+(1U<<0));// C is the irreducible poly. (64,4,3,1,0)
    __m128i Q2 = _mm_clmulepi64_si128( A, C, 0x01);
    __m128i Q3 = _mm_shuffle_epi8(_mm_setr_epi8(0,  27,  54,  45,  108,  119,  90,  65,  216,  195,  238,  245,  180,  175,  130,  153),
    				_mm_srli_si128(Q2,8));
    __m128i Q4 =  _mm_xor_si128(Q2,A);
    const __m128i final  = _mm_xor_si128(Q3,Q4);
    return final;/// WARNING: HIGH 64 BITS CONTAIN GARBAGE
 }



uint64_t precompReduction64( __m128i A) {
	return _mm_cvtsi128_si64(precompReduction64_si128(A));
}







__m128i barrettWithoutPrecomputation16_si128( __m128i A) {
    ///http://www.jjj.de/mathdata/minweight-primpoly.txt
    const uint64_t irredpoly = 1UL+(1UL<<2)+(1UL<<3)+(1UL<<5)+(1UL<<16);
    // it is important, for the algo. we have chosen that 5 is smaller
    // equal than 8=16/2
    const int n = 16;// degree of the polynomial
    const __m128i C = _mm_set_epi64x(0,irredpoly);// C is the irreducible poly.
    const __m128i Q1 = _mm_srli_si128 (A, 2);
    const __m128i Q2 = _mm_clmulepi64_si128( Q1, C, 0x00);// A div x^n
    const __m128i Q3 = _mm_srli_si128 (Q2, 2);
    const __m128i Q4 = _mm_clmulepi64_si128( Q3, C, 0x00);
    const __m128i final  = _mm_xor_si128 (A, Q4);
    return final;
}

uint16_t barrettWithoutPrecomputation16( __m128i A) {
    return (uint16_t) _mm_cvtsi128_si32(barrettWithoutPrecomputation16_si128(A));
}









#endif
#endif
