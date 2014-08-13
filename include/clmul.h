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




// multiplication with lazy reduction
// assumes that the two highest bits of one of the inputs are zero
// returns a lazy reduction
__m128i mul128by128to128_lazymod127( __m128i A, __m128i B) {
	__m128i Alow = _mm_clmulepi64_si128(A,B,0x00);
	__m128i Ahigh = _mm_clmulepi64_si128(A,B,0x11);
	__m128i Amix1 = _mm_clmulepi64_si128(A,B,0x01);
	__m128i Amix2 = _mm_clmulepi64_si128(A,B,0x10);
	__m128i Amix = _mm_xor_si128(Amix1,Amix2);
	Ahigh = _mm_xor_si128(Ahigh,_mm_srli_si128(Amix,8));
	Alow = _mm_xor_si128(Alow,_mm_slli_si128(Amix,8));
	// now the lazy reduction
	///////////////////////////////////////////////////
	// We want to take Ahigh and compute       (  Ahigh <<1 ) XOR (  Ahigh <<2 )
	// Except that there is no way to shift an entire XMM register by 1 or 2 bits  using a single instruction.
	// So how do you compute Ahigh <<1 using as few instructions as possible?
	//
	// First you do _mm_slli_epi64(Ahigh,1). This is *almost* correct... except that
	// the 64th bit is not shifted in 65th position.
	// Well, ok, but we can compute Ahigh >> 8, this given by _mm_srli_si128(Ahigh,1)
	// _mm_srli_si128(Ahigh,1) has the same bits as Ahigh (except that we lose the lowest 8)
	// but at different positions.
	// So let us shift left the result again...
	//  _mm_slli_epi64(_mm_srli_si128(Ahigh,1),1)
	// If you keep track, this is "almost" equivalent to A >> 7, except that the 72th bit
	// from A is lost.
	// From something that is almost A >>7, we can get back something that is almost A << 1
	// by shifting left by 8 bits...
	// _mm_slli_si128(_mm_slli_epi64(_mm_srli_si128(Ahigh,1),1),1)
	// So this is a shift left by 1, except that the 72th bit is lost along with the lowest 8 bits.
	// We have that  _mm_slli_epi64(Ahigh,1) is a shift let by 1 except that the 64th bit
	// is lost. We can combine the two to get the desired result (just OR them).
	// The algorithm below is just an optimized version of this where we do both shifts (by 1 and 2)
	// at the same time and XOR the result.
	//
	__m128i shifteddownAhigh = _mm_srli_si128(Ahigh,1);
	__m128i s1 = _mm_slli_epi64(Ahigh,1);
	__m128i s2 = _mm_slli_epi64(Ahigh,2);
	__m128i sd1 = _mm_slli_si128(_mm_slli_epi64(shifteddownAhigh,1),1);
	__m128i sd2 = _mm_slli_si128(_mm_slli_epi64(shifteddownAhigh,2),1);
	s1 = _mm_or_si128(s1,sd1);
	s2 = _mm_or_si128(s2,sd2);
	__m128i reduced = _mm_xor_si128(s1,s2);
	// combining results
	__m128i final = _mm_xor_si128(Alow,reduced);
	return reduced;
}



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
	const __m128i final = _mm_xor_si128 (A, Q4);
	return final;/// WARNING: HIGH 96 BITS CONTAIN GARBAGE
}

uint32_t barrettWithoutPrecomputation32( __m128i A) {
	return _mm_cvtsi128_si32(barrettWithoutPrecomputation32_si128(A));
}


// modulo reduction to 64-bit value. The high 64 bits contain garbage, see precompReduction64
__m128i precompReduction64_si128( __m128i A) {

	const __m128i C = _mm_set_epi64x(1U,(1U<<4)+(1U<<3)+(1U<<1)+(1U<<0)); // C is the irreducible poly. (64,4,3,1,0)
	__m128i Q2 = _mm_clmulepi64_si128( A, C, 0x01);
	__m128i Q3 = _mm_shuffle_epi8(_mm_setr_epi8(0, 27, 54, 45, 108, 119, 90, 65, 216, 195, 238, 245, 180, 175, 130, 153),
			_mm_srli_si128(Q2,8));
	__m128i Q4 = _mm_xor_si128(Q2,A);
	const __m128i final = _mm_xor_si128(Q3,Q4);
	return final;/// WARNING: HIGH 64 BITS CONTAIN GARBAGE
}

uint64_t precompReduction64( __m128i A) {
	return _mm_cvtsi128_si64(precompReduction64_si128(A));
}



#endif
#endif
