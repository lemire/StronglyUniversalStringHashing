#ifndef _CLMUL_H_
#define _CLMUL_H_
#ifdef __AVX__ // intel does not define PCLMUL
#define __PCLMUL__ 1
#endif

#ifdef __PCLMUL__

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

//////////////////////////
/// next line would be the "right thing to do" but instead we include
/// just the SSE/AVX dependencies we need otherwise the fact that
/// we mix C and C++ leads to a build error with Intel compilers.
/// See https://github.com/lemire/StronglyUniversalStringHashing/issues/19
//////////////////////////
//#include <x86intrin.h> // can't do that on Intel compiler.
#include <wmmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
///////// End of compatibility hack.

void printme32(__m128i v1) {
    printf(" %u %u %u %u  ", _mm_extract_epi32(v1,0), _mm_extract_epi32(v1,1), _mm_extract_epi32(v1,2), _mm_extract_epi32(v1,3));
}

void printme64(__m128i v1) {
    printf(" %llu %llu  ", _mm_extract_epi64(v1,0), _mm_extract_epi64(v1,1));
}

// computes a << 1
static inline __m128i leftshift1(__m128i a) {
    const int x = 1;
    __m128i u64shift =  _mm_slli_epi64(a,x);
    __m128i topbits =  _mm_slli_si128(_mm_srli_epi64(a,64 - x),sizeof(uint64_t));
    return _mm_or_si128(u64shift, topbits);
}

// computes a << 2
static inline __m128i leftshift2(__m128i a) {
    const int x = 2;
    __m128i u64shift =  _mm_slli_epi64(a,x);
    __m128i topbits =  _mm_slli_si128(_mm_srli_epi64(a,64 - x),sizeof(uint64_t));
    return _mm_or_si128(u64shift, topbits);
}

//////////////////
// compute the "lazy" modulo with 2^127 + 2 + 1, actually we compute the
// modulo with (2^128 + 4 + 2) = 2 * (2^127 + 2 + 1) ,
// though  (2^128 + 4 + 2) is not
// irreducible, we have that
//     (x mod (2^128 + 4 + 2)) mod (2^127 + 2 + 1) == x mod (2^127 + 2 + 1)
// That's true because, in general ( x mod k y ) mod y = x mod y.
//
// Precondition:  given that Ahigh|Alow represents a 254-bit value
//                  (two highest bits of Ahigh must be zero)
//////////////////
__m128i lazymod127(__m128i Alow, __m128i Ahigh) {
    ///////////////////////////////////////////////////
    // CHECKING THE PRECONDITION:
    // Important: we are assuming that the two highest bits of Ahigh
    // are zero. This could be checked by adding a line such as this one:
    // if(_mm_extract_epi64(Ahigh,1) >= (1ULL<<62)){printf("bug\n");abort();}
    //                       (this assumes SSE4.1 support)
    ///////////////////////////////////////////////////
    // The answer is Alow XOR  (  Ahigh <<1 ) XOR (  Ahigh <<2 )
    // This is correct because the two highest bits of Ahigh are
    // assumed to be zero.
    ///////////////////////////////////////////////////
    // credit for simplified implementation : Jan Wassenberg
    __m128i shift1 = leftshift1(Ahigh);
    __m128i shift2 = leftshift2(Ahigh);
    __m128i final =  _mm_xor_si128(_mm_xor_si128(Alow, shift1),shift2);
    return final;
}


// multiplication with lazy reduction
// assumes that the two highest bits of the 256-bit multiplication are zeros
// returns a lazy reduction
__m128i mul128by128to128_lazymod127( __m128i A, __m128i B) {
    __m128i Amix1 = _mm_clmulepi64_si128(A,B,0x01);
    __m128i Amix2 = _mm_clmulepi64_si128(A,B,0x10);
    __m128i Alow = _mm_clmulepi64_si128(A,B,0x00);
    __m128i Ahigh = _mm_clmulepi64_si128(A,B,0x11);
    __m128i Amix = _mm_xor_si128(Amix1,Amix2);
    Amix1 = _mm_slli_si128(Amix,8);
    Amix2 = _mm_srli_si128(Amix,8);
    Alow = _mm_xor_si128(Alow,Amix1);
    Ahigh = _mm_xor_si128(Ahigh,Amix2);
    return lazymod127(Alow, Ahigh);
}

// multiplication with lazy reduction
// A1 * B1 + A2 * B2
// assumes that the two highest bits of the 256-bit multiplication are zeros
// returns a lazy reduction
__m128i mul128by128to128_lazymod127_2by2( __m128i A1,__m128i A2,
        __m128i B1,__m128i B2) {
    __m128i Amix11 = _mm_clmulepi64_si128(A1,B1,0x01);
    __m128i Amix21 = _mm_clmulepi64_si128(A1,B1,0x10);
    __m128i Amix12 = _mm_clmulepi64_si128(A2,B2,0x01);
    __m128i Amix22 = _mm_clmulepi64_si128(A2,B2,0x10);
    __m128i Alow1 = _mm_clmulepi64_si128(A1,B1,0x00);
    __m128i Ahigh1 = _mm_clmulepi64_si128(A1,B1,0x11);
    __m128i Alow2 = _mm_clmulepi64_si128(A2,B2,0x00);
    __m128i Ahigh2 = _mm_clmulepi64_si128(A2,B2,0x11);

    __m128i Amix1 = _mm_xor_si128(Amix11,Amix21);
    __m128i Amix2 = _mm_xor_si128(Amix12,Amix22);
    __m128i Amix = _mm_xor_si128(Amix1,Amix2);
    Amix1 = _mm_slli_si128(Amix,8);
    Amix2 = _mm_srli_si128(Amix,8);
    __m128i Alow = _mm_xor_si128(Alow1,Alow2);
    __m128i Ahigh = _mm_xor_si128(Ahigh1,Ahigh2);
    Alow = _mm_xor_si128(Alow,Amix1);
    Ahigh = _mm_xor_si128(Ahigh,Amix2);
    return lazymod127(Alow, Ahigh);
}


// multiplication with lazy reduction
// A1 * B1 + A2 * B2 +  A3 * B3 + A4 * B4
// assumes that the two highest bits of the 256-bit multiplication are zeros
// returns a lazy reduction
__m128i mul128by128to128_lazymod127_4by4( __m128i A1,__m128i A2,__m128i A3,__m128i A4,
        __m128i B1,__m128i B2,__m128i B3,__m128i B4) {
    __m128i Amix11 = _mm_clmulepi64_si128(A1,B1,0x01);
    __m128i Amix21 = _mm_clmulepi64_si128(A1,B1,0x10);
    __m128i Amix12 = _mm_clmulepi64_si128(A2,B2,0x01);
    __m128i Amix22 = _mm_clmulepi64_si128(A2,B2,0x10);
    __m128i Amix13 = _mm_clmulepi64_si128(A3,B3,0x01);
    __m128i Amix23 = _mm_clmulepi64_si128(A3,B3,0x10);
    __m128i Amix14 = _mm_clmulepi64_si128(A4,B4,0x01);
    __m128i Amix24 = _mm_clmulepi64_si128(A4,B4,0x10);

    __m128i Alow1 = _mm_clmulepi64_si128(A1,B1,0x00);
    __m128i Ahigh1 = _mm_clmulepi64_si128(A1,B1,0x11);
    __m128i Alow2 = _mm_clmulepi64_si128(A2,B2,0x00);
    __m128i Ahigh2 = _mm_clmulepi64_si128(A2,B2,0x11);
    __m128i Alow3 = _mm_clmulepi64_si128(A3,B3,0x00);
    __m128i Ahigh3 = _mm_clmulepi64_si128(A3,B3,0x11);
    __m128i Alow4 = _mm_clmulepi64_si128(A4,B4,0x00);
    __m128i Ahigh4 = _mm_clmulepi64_si128(A4,B4,0x11);



    __m128i Amix1 = _mm_xor_si128(Amix11,Amix21);
    __m128i Amix2 = _mm_xor_si128(Amix12,Amix22);
    __m128i Amix3 = _mm_xor_si128(Amix13,Amix23);
    __m128i Amix4 = _mm_xor_si128(Amix14,Amix24);

    Amix12 = _mm_xor_si128(Amix1,Amix2);
    Amix23 = _mm_xor_si128(Amix3,Amix4);
    __m128i Amix = _mm_xor_si128(Amix12,Amix23);

    Amix1 = _mm_slli_si128(Amix,8);

    Amix2 = _mm_srli_si128(Amix,8);

    __m128i Alow12 = _mm_xor_si128(Alow1,Alow2);
    __m128i Alow34 = _mm_xor_si128(Alow3,Alow4);
    __m128i Alow = _mm_xor_si128(Alow12,Alow34);

    __m128i Ahigh12 = _mm_xor_si128(Ahigh1,Ahigh2);
    __m128i Ahigh34 = _mm_xor_si128(Ahigh3,Ahigh4);
    __m128i Ahigh = _mm_xor_si128(Ahigh12,Ahigh34);

    Alow = _mm_xor_si128(Alow,Amix1);
    Ahigh = _mm_xor_si128(Ahigh,Amix2);
    return lazymod127(Alow, Ahigh);
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
    //const int n = 32;// degree of the polynomial
    const __m128i C = _mm_cvtsi64_si128(irredpoly);// C is the irreducible poly.
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

    //const __m128i C = _mm_set_epi64x(1U,(1U<<4)+(1U<<3)+(1U<<1)+(1U<<0)); // C is the irreducible poly. (64,4,3,1,0)
    const __m128i C = _mm_cvtsi64_si128((1U<<4)+(1U<<3)+(1U<<1)+(1U<<0));
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
