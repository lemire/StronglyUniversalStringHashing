#include <stdio.h>
#include "clmul.h"

#ifdef __PCLMUL__

int is_zero(__m128i a){
	return _mm_testz_si128(a,a);
}

// checks to see if a == b
int equal(__m128i a, __m128i b) {
	__m128i xorm = _mm_xor_si128(a,b);
	return _mm_testz_si128(xorm,xorm);

}

void printme32(__m128i v1) {
    printf(" %u %u %u %u  ", _mm_extract_epi32(v1,0), _mm_extract_epi32(v1,1), _mm_extract_epi32(v1,2), _mm_extract_epi32(v1,3));
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
		a = _mm_xor_si128(a,rmult);
		int newdegreea = degree(a);
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



void displayfirst() {
	for(uint64_t a = 0; a< 16; ++a) {
		const __m128i A = _mm_set_epi64x(a,0);
        __m128i fastmod = barrettWithoutPrecomputation64_si128(A);
        printf(" %llu, ",_mm_extract_epi64(fastmod,0));

	}
	printf("\n");

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

void precompclmulunittest0_64() {
    printf("CLMUL test precomp 0_64...\n");
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
        __m128i fastmod = precompReduction64_si128(A);
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
	uint64_t mul2 = 12344567788889+(1ULL<<62)+(1ULL<<61)+(1ULL<<55);//random-like
    for(uint64_t a = 1; a< 1024; ++a) {
        const __m128i A = _mm_set_epi64x(0,mul2*a);
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
void clmulunittest0_16() {
    printf("CLMUL test 0_16...\n");
    const uint64_t irredpoly = 1UL+(1UL<<2)+(1UL<<3)+(1UL<<5)+(1UL<<16);
    const __m128i C = _mm_set_epi64x(0,irredpoly);// C is the irreducible poly.
	uint32_t mul2 = 2979263833U+(1U<<31);//random-like
    for(uint64_t a = 1; a< 1024; ++a) {
        const __m128i A = _mm_set_epi32(0,0,0,mul2*a);
        __m128i sillymod = slowcarrylessdivision(A,A).y;
        if(!is_zero(sillymod)) {
            printme64(sillymod);
        	printf("silly mod is not zero?\n");
        	abort();
        }
        __m128i slowmod = slowcarrylessdivision(A,C).y;
        __m128i fastmod = barrettWithoutPrecomputation16_si128(A);
        fastmod = _mm_and_si128(fastmod,_mm_set_epi32(0,0,0,0xFFFF));// keep just the low 16 bits
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
    //displayfirst();
    precompclmulunittest0_64();
    clmulunittest0_64();
    clmulunittest0_32();
    clmulunittest0_16();
    clmulunittest1();
    clmulunittest2();
    clmulunittest3();
    clmulunittest3a();
    printf("CLMUL code looks ok.\n");
}


int main() {
	clmulunittests();
	return 0;
}

#else
int main() {
	printf("clmul instruction set not detected. Testing not done.");
	return 0;
}

#endif

