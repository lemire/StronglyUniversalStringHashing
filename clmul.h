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



int is_zero(__m128i a){
	return _mm_testz_si128(a,a);
}

// checks to see if a == b
int equal(__m128i a, __m128i b) {
	__m128i xorm = _mm_xor_si128(a,b);
	return _mm_testz_si128(xorm,xorm);

}

void printme32(__m128i v1) {
    printf(" %lu %lu %lu %lu  ", _mm_extract_epi32(v1,0), _mm_extract_epi32(v1,1), _mm_extract_epi32(v1,2), _mm_extract_epi32(v1,3));
}

void printme64(__m128i v1) {
    printf(" %llu %llu  ", _mm_extract_epi64(v1,0), _mm_extract_epi64(v1,1));
}

// useful for debugging, essentially determines the most significant bit set
// return 0 if a is made solely of zeroes
int degree(__m128i a) {
	int x4 = _mm_extract_epi32 (a,3);
	if(x4 != 0) {
		return  32 - __builtin_clz(x4) + 3 * 32 - 1;
	}
	int x3 = _mm_extract_epi32 (a,2);
	if(x3 != 0) {
		return 32 - __builtin_clz(x3) + 2 * 32 - 1;
	}
	int x2 = _mm_extract_epi32 (a,1);
	if(x2 != 0) {
		return 32 - __builtin_clz(x2) +  32 - 1;
	}
	int x1 = _mm_extract_epi32 (a,0);
	if(x1 != 0) {
		return 32 - __builtin_clz(x1) - 1;
	}
	return 0;
}

typedef struct  {
	__m128i x;
	__m128i y;
} PairOfVec;

// 128bit x 128bit to 128bit multiplication
__m128i fullcarrylessmultiply(__m128i a, __m128i b) {
	__m128i part1 = _mm_clmulepi64_si128( a, b, 0x00);
	__m128i part2 = _mm_clmulepi64_si128( a, b, 0x10);
	part2 = _mm_slli_si128(part2,8);
	__m128i part3 = _mm_clmulepi64_si128( a, b, 0x01);
	part3 = _mm_slli_si128(part3,8);
	return _mm_xor_si128( _mm_xor_si128(part1,part2),part3);

}
// useful for debugging. Divides a by b in carryless model
// surely, there is a better way to do this...
PairOfVec slowcarrylessdivision(__m128i a, __m128i b) {
	__m128i copya = a;
	int degreea = degree(a);
	int degreeb = degree(b);
	PairOfVec result ;
	result.x = _mm_setzero_si128();
	printf("degreea = %i degreeb = %i \n",degreea,degreeb);
	while(degreea > 0) {
		int off = degreea - degreeb;
		if(off<0) break;
		__m128i bitf = off>=64 ? _mm_set_epi64x(1ULL<<(off-64),0) : _mm_set_epi64x(0,1ULL<<off);
		int degreebitf = degree(bitf);
		if( degreebitf + degreeb != degreea) {
			printme64(bitf);
			printf("you would think that this would work %i %i %i %i \n", degreebitf,degreeb,degreea,off);
		}
		result.x = _mm_xor_si128(bitf,result.x);
		__m128i rmult = fullcarrylessmultiply(b, bitf);
		//printf("result of the multiplication %i %i  %i %i \n",degree(rmult),degreea,degreeb,degreebitf);
		a = _mm_xor_si128(a,rmult);
		int newdegreea = degree(a);
		//printf("newdegreea = %i \n",newdegreea);
		// printme64(a);
		if(newdegreea >= degreea)
			printf("this will not end well\n");
		degreea = newdegreea;
	}
	result.y = a;
	// now we check our result
	if(!equal(_mm_xor_si128(fullcarrylessmultiply(result.x, b),result.y),copya)){
		printf("Slow division is buggy\n");
		abort();

	};

	return result;
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
    const __m128i Q1 = _mm_srli_epi64 (A, n);
    const __m128i Q2 = _mm_clmulepi64_si128( Q1, C, 0x00);// A div x^n
    const __m128i Q3 = _mm_srli_epi64 (Q2, n);
    // commenting out the long way derived from the paper (following two lines are enough)
    //__m128i R1 = _mm_and_si128 (maskm128,A);
    //__m128i R2 = _mm_and_si128 (maskm128,_mm_clmulepi64_si128( Q3, C, 0x00));
    //__m128i final  = _mm_xor_si128 (R1, R2);
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


// standard unoptimized 32-bit CLMUL hashing
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


// optimized 32-bit CLMUL hashing
uint32_t hashGaloisFieldMultilinearHalfMultiplications(const uint64_t*  randomsource, const uint32_t *  string, const size_t length) {
    const uint32_t * const endstring = string + length;
    const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    __m128i acc = _mm_set_epi64x(0,*(randomsource32));
    randomsource32 += 1;
    for(; string +3 < endstring; randomsource32+=4,string+=4 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource32);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i twosums = _mm_xor_si128(temp1,temp2); 
        const __m128i part1 = _mm_unpacklo_epi32(twosums,_mm_setzero_si128());
        const __m128i clprod1  = _mm_clmulepi64_si128( part1, part1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);   
        const __m128i part2 = _mm_unpackhi_epi32(twosums,_mm_setzero_si128());
        const __m128i clprod2  = _mm_clmulepi64_si128( part2, part2, 0x10);
        acc = _mm_xor_si128 (clprod2,acc);   
     }
    if(string + 1  < endstring) {
        __m128i temp1 = _mm_set_epi64x(*randomsource32,*(randomsource32+1));
        __m128i temp2 = _mm_set_epi64x(*string,*(string+1));
        __m128i twosums = _mm_xor_si128(temp1,temp2);
        __m128i clprod  = _mm_clmulepi64_si128( twosums, twosums, 0x10);
        acc = _mm_xor_si128 (clprod,acc);
        randomsource32+=2;
        string+=2;
    }
    if( string!= endstring ) {
        __m128i temp = _mm_set_epi64x(*randomsource32,*string);
        __m128i clprod  = _mm_clmulepi64_si128( temp, temp, 0x10);
        acc = _mm_xor_si128 (clprod,acc);
    }
    return barrettWithoutPrecomputation32(acc);
}


// a 64-bit version
// strongly universal and regular
uint64_t hashGaloisFieldfast64(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
    const uint64_t * const endstring = string + length;
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
uint64_t hashGaloisFieldfast64half(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
    const uint64_t * const endstring = string + length;
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
uint64_t hashGaloisFieldfast64halfunrolled(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length;
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
uint64_t hashGaloisFieldPoly64(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
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


// fast 64-bit polynomial hashing, uses only one key
// expected to be fast!
//TODO: can use more keys for increased universality
uint64_t fasthashGaloisFieldPoly64(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
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
        IACA_START;
        __m128i temp = _mm_lddqu_si128((__m128i *) string);
        __m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2));
        const __m128i x1 =  _mm_and_si128 (temp,mask);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp, key, 0x01);
        const __m128i clprod2  = _mm_clmulepi64_si128( temp2, key, 0x10);
        const __m128i clprod3  = _mm_clmulepi64_si128( temp2, key2, 0x01);
        acc  = _mm_clmulepi64_si128( acc, key2, 0x10);        
        acc = _mm_xor_si128 (acc,_mm_xor_si128 (_mm_xor_si128 (x1,clprod1),_mm_xor_si128 (clprod2,clprod3)));
        IACA_END;
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

// like MHH, this is essentially multilinear with 64bit multiplication
// summed up over a 128-bit counter 
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





void clmulunittest0_64() {
    printf("CLMUL test 0_64...\n");
    const __m128i C = _mm_set_epi64x(1U,(1U<<4)+(1U<<3)+(1U<<1)+(1U<<0));// C is the irreducible poly. (64,4,3,1,0)
	uint64_t mul1 = 4343232+(1ULL<<63)+(1ULL<<60)+(1ULL<<45);//random-like
	uint64_t mul2 = 12344567788889+(1ULL<<62)+(1ULL<<61)+(1ULL<<55);//random-like
    for(uint64_t a = 1; a< 1024; ++a) {
        const __m128i A = _mm_set_epi64x(mul1*a,mul2*a);
        __m128i sillymod = slowcarrylessdivision(A,A).y;
        if(!is_zero(sillymod)) {
        	  printme64(sillymod);
              printf("silly mod is not zero?\n");
              abort();
        }
        PairOfVec AdivC = slowcarrylessdivision(A,C);
        __m128i slowmod = AdivC.y;
        __m128i fastmod = barrettWithoutPrecomputation64_si128(A);
        fastmod = _mm_and_si128(fastmod,_mm_set_epi64x(0,-1));// keep just the low 64 bits

        if(!equal(slowmod,fastmod)) {
        	printf("slowmod = ");
        	printme64(slowmod);
      	    printf("\n");
      	    printf("fastmod = ");
    	    printme64(fastmod);
    	    printf("\n");
    	    printf("64-bit bug slowmod and fastmod differs\n");
    	    abort();
        }
    }
    printf("Test passed!\n");

}


void clmulunittest0_32() {
    printf("CLMUL test 0_32...\n");
    const uint64_t irredpoly = 1UL+(1UL<<2)+(1UL<<6)+(1UL<<7)+(1UL<<32);
    const __m128i C = _mm_set_epi64x(0,irredpoly);// C is the irreducible poly.
	uint64_t mul1 = 4343232+(1ULL<<63)+(1ULL<<60)+(1ULL<<45);//random-like
	uint64_t mul2 = 12344567788889+(1ULL<<62)+(1ULL<<61)+(1ULL<<55);//random-like
    for(uint64_t a = 1; a< 1024; ++a) {
        const __m128i A = _mm_set_epi64x(mul1*a,mul2*a);
        __m128i sillymod = slowcarrylessdivision(A,A).y;
        if(!is_zero(sillymod)) {
            printme64(sillymod);
        	printf("silly mod is not zero?\n");
        	abort();
        }
        __m128i slowmod = slowcarrylessdivision(A,C).y;
        __m128i fastmod = barrettWithoutPrecomputation32_si128(A);
        fastmod = _mm_and_si128(fastmod,_mm_set_epi32(0,0,0,-1));// keep just the low 32 bits

        if(!equal(slowmod,fastmod)) {
        	printf("slowmod = ");
        	printme32(slowmod);
        	printf("\n");
        	printf("fastmod = ");
        	printme32(fastmod);
        	printf("\n");
        	printf("32-bit bug slowmod and fastmod differs\n");
        	abort();
        }
    }
    printf("Test passed!\n");

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
    clmulunittest0_64();
    clmulunittest0_32();
    clmulunittest1();
    clmulunittest2();
    clmulunittest3();
    clmulunittest3a();
    printf("CLMUL code looks ok.\n");
}

#endif

#endif
