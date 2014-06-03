#ifndef _CLMUL_H_
#define _CLMUL_H_


#ifdef __PCLMUL__

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <x86intrin.h>


#define IACA
#ifdef IACA
#include </opt/intel/iaca/include/iacaMarks.h>
#endif


void printme(__m128i v1) {
    printf(" %i %i %i %i times ", _mm_extract_epi32(v1,0), _mm_extract_epi32(v1,1), _mm_extract_epi32(v1,2), _mm_extract_epi32(v1,3));
}

void printme64(__m128i v1) {
    printf(" %llu %llu times ", _mm_extract_epi64(v1,0), _mm_extract_epi64(v1,1));
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
uint32_t barrettWithoutPrecomputation32( __m128i A) {
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
    const __m128i Q1 = _mm_srli_epi64 (A, n);
    const __m128i Q2 = _mm_clmulepi64_si128( Q1, C, 0x00);// A div x^n
    const __m128i Q3 = _mm_srli_epi64 (Q2, n);
    // commenting out the long way derived from the paper (following two lines are enough)
    //__m128i R1 = _mm_and_si128 (maskm128,A);
    //__m128i R2 = _mm_and_si128 (maskm128,_mm_clmulepi64_si128( Q3, C, 0x00));
    //__m128i final  = _mm_xor_si128 (R1, R2);
    const __m128i Q4 = _mm_clmulepi64_si128( Q3, C, 0x00);
    const __m128i final  = _mm_xor_si128 (A, Q4);
    return _mm_cvtsi128_si32(final);
}
uint64_t barrettWithoutPrecomputation64( __m128i A) {
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
    //printf ("A = "); printme64(A); printf("\n");

    const __m128i Q1 = _mm_srli_si128 (A, 8);
    //printf ("Q1 = "); printme64(Q1); printf("\n");
    __m128i Q2 = _mm_clmulepi64_si128( Q1, C, 0x00);
    Q2 = _mm_xor_si128 (Q2,_mm_slli_si128 (Q1, 8) );
    //printf ("Q2 = "); printme64(Q2); printf("\n");
    const __m128i Q3 = _mm_srli_si128 (Q2, 8);
    const __m128i Q4 = _mm_clmulepi64_si128( Q3, C, 0x00);
    const __m128i final  = _mm_xor_si128 (A, Q4);
    //printf ("final = "); printme64(final); printf("\n");


    return _mm_cvtsi128_si64(final);
}


uint16_t barrettWithoutPrecomputation16( __m128i A) {
    ///http://www.jjj.de/mathdata/minweight-primpoly.txt
    const uint64_t irredpoly = 1UL+(1UL<<2)+(1UL<<3)+(1UL<<5)+(1UL<<16);
    // it is important, for the algo. we have chosen that 5 is smaller
    // equal than 8=16/2
    const int n = 16;// degree of the polynomial
    const __m128i C = _mm_set_epi64x(0,irredpoly);// C is the irreducible poly.
    const __m128i Q1 = _mm_srli_epi64 (A, n);
    const __m128i Q2 = _mm_clmulepi64_si128( Q1, C, 0x00);// A div x^n
    const __m128i Q3 = _mm_srli_epi64 (Q2, n);
    const __m128i Q4 = _mm_clmulepi64_si128( Q3, C, 0x00);
    const __m128i final  = _mm_xor_si128 (A, Q4);
    return (uint16_t) _mm_cvtsi128_si32(final);
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
    return barrettWithoutPrecomputation32(acc);
}



uint32_t hashGaloisFieldMultilinearHalfMultiplications(const uint64_t*  randomsource, const uint32_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint32_t * const endstring = string + length;
    const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    __m128i acc = _mm_set_epi64x(0,*(randomsource32));
    randomsource32 += 1;
    for(; string!= endstring; randomsource32+=2,string+=2 ) {
#ifdef IACA
        IACA_START;// place after for loop
#endif
        __m128i temp1 = _mm_set_epi64x(*randomsource32,*(randomsource32+1));
        __m128i temp2 = _mm_set_epi64x(*string,*(string+1));
        __m128i twosums = _mm_xor_si128(temp1,temp2);
        __m128i clprod  = _mm_clmulepi64_si128( twosums, twosums, 0x10);
        acc = _mm_xor_si128 (clprod,acc);
#ifdef IACA
        IACA_END;// place after loop
#endif
    }
    return barrettWithoutPrecomputation32(acc);
}


// a 64-bit version
uint64_t hashGaloisFieldfast64(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length/2*2;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        const __m128i clprod2  = _mm_clmulepi64_si128( temp1, temp2, 0x11);
        acc = _mm_xor_si128 (clprod1,acc);
        acc = _mm_xor_si128 (clprod2,acc);
    }
    assert(string == endstring);
    return barrettWithoutPrecomputation64(acc);
}


// a 64-bit version with half the number of multiplications
uint64_t hashGaloisFieldfast64half(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length*2/2;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    assert(string == endstring);
    return barrettWithoutPrecomputation64(acc);
}


// a 64-bit version with half the number of multiplications
uint64_t hashGaloisFieldfast64halfunrolled(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length*2/2;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string!= endstring; randomsource+=4,string+=4 ) {
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
    assert(string == endstring);
    return barrettWithoutPrecomputation64(acc);
}


// like MHH
uint64_t referenceproduct(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {

    uint64_t low = 0;
    uint64_t high = 0;
    size_t i = 0;
    for(; i<length/8*8; i+=8) {
    __asm__ (
    "movq (%[u]),%%rax\n"
    "mulq (%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 8(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 8(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 16(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 16(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 24(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 24(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 32(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 32(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 40(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 40(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 48(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 48(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 56(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 56(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "adcq %%rdx,  %[rh]\n"
                 :  [rh] "+r" (high) , [rl] "+r" (low)  : [u] "r" (randomsource+i), [v] "r" (string+i)  :"rdx","rax","memory","cc");
	}

    for(; i<length; ++i) {
        __asm__ ("mulq %[v]\n"
                 "addq %%rax,  %[rl]\n"
                 "adcq %%rdx,  %[rh]\n"
                 :  [rh] "+r" (high), [rl] "+r" (low)  : [u] "a" (randomsource[i]), [v] "r" (string[i])  :"rdx","cc");
    }
    return low+high;// should be modulo
}








void clmulunittest1() {
    printf("CLMUL test 1...\n");
    // A * x = y ought to be invertible.
    // brute force check is hard, but we can do some checking
    // in the 16-bit case
    for(int a = 1; a< 1<<16; a+=32) {
        __m128i A = _mm_set_epi32(0,0,0,2*a);
        int counter[1<<16];
        for(int j = 0; j < 1<<16; ++j) counter[j] = 0;
        for(int k = 0; k< 1<<16; ++k) {
            __m128i v1 = _mm_set_epi16(0,0,0,0,0,0,0,k);
            __m128i clprod1  = _mm_clmulepi64_si128( A, v1, 0x00);
            uint64_t m64 = barrettWithoutPrecomputation64(clprod1);
            uint32_t m32 = barrettWithoutPrecomputation32(_mm_set_epi64x(0,m64));
            uint16_t m16 = barrettWithoutPrecomputation16(_mm_set_epi32(0,0,0,m32));
            counter[m16]++;
        }
        for(int j = 0; j < 1<<16; ++j) {
            if(counter[j] != 1) {
                printf("j= %i c = %i\n",j,counter[j]);
                printf("bug\n");
                abort();
            }

        }
    }
}




void clmulunittest2() {
    printf("CLMUL test 2...\n");
    // Idea: we fish for two non-zero 32-bit integers with a zero product
    for(int a = 1; a< 1<<16; a+=32) {
        __m128i A = _mm_set_epi32(0,0,0,2*a);
        for(int k = 1; k< 1<<16; ++k) {
            __m128i v1 = _mm_set_epi16(0,0,0,0,0,0,k,0);
            __m128i clprod1  = _mm_clmulepi64_si128( A, v1, 0x00);
            uint64_t m64 = barrettWithoutPrecomputation64(clprod1);
            uint32_t m32 = barrettWithoutPrecomputation32(_mm_set_epi64x(0,m64));
            if(m32 == 0) {
                printf("bug\n");
                abort();

            }
        }
    }
}

void clmulunittest3() {
    printf("CLMUL test 3...\n");
    // Idea: we fish for two non-zero 64-bit integers with a zero product
    for(int a = 1; a< 1<<16; a+=32) {
        __m128i A = _mm_set_epi32(0,0,0,2*a);
        for(int k = 1; k< 1<<16; ++k) {
            __m128i v1 = _mm_set_epi16(0,0,0,0,k,0,0,0);
            __m128i clprod1  = _mm_clmulepi64_si128( A, v1, 0x00);
            uint64_t m64 = barrettWithoutPrecomputation64(clprod1);
            if(m64 == 0) {
                printf("A=");
                printme64(v1);
                printf("\n");
                printf("v1=");
                printme64(v1);
                printf("\n");
                printf("clprod1=");
                printme64(clprod1);
                printf("\n");
                printf("bug\n");
                abort();
            }
        }
    }
}
void clmulunittest3a() {
    printf("CLMUL test 3a...\n");
    // Idea: we fish for two non-zero 64-bit integers with a zero product
    for(int a = 1; a< 1<<16; a+=32) {
        __m128i A = _mm_set_epi32(0,0,2*a,0);
        for(int k = 1; k< 1<<16; ++k) {
            __m128i v1 = _mm_set_epi16(0,0,0,0,k,0,0,0);
            __m128i clprod1  = _mm_clmulepi64_si128( A, v1, 0x00);
            uint64_t m64 = barrettWithoutPrecomputation64(clprod1);
            if(m64 == 0) {
                printf("A=");
                printme64(v1);
                printf("\n");
                printf("v1=");
                printme64(v1);
                printf("\n");
                printf("clprod1=");
                printme64(clprod1);
                printf("\n");
                printf("bug\n");
                abort();
            }
        }
    }
}


void clmulunittests() {
    printf("Testing CLMUL code...\n");
    clmulunittest1();
    clmulunittest2();
    clmulunittest3();
    clmulunittest3a();
    printf("CLMUL code looks ok.\n");
}

#endif

#endif
