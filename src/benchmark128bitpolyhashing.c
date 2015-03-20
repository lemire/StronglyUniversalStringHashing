/*
 * Benchmarking 128-bit polynomial hashing
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>

#ifdef __AVX__
#define __PCLMUL__ 1
#endif
#include "clmul.h"


typedef unsigned long long ticks;

// Taken from stackoverflow (see http://stackoverflow.com/questions/3830883/cpu-cycle-count-based-profiling-in-c-c-linux-x86-64)
// Can give nonsensical results on multi-core AMD processors.
ticks rdtsc() {
	unsigned int lo, hi;
	asm volatile (
			"cpuid \n" /* serializing */
			"rdtsc"
			: "=a"(lo), "=d"(hi) /* outputs */
			: "a"(0) /* inputs */
			: "%ebx", "%ecx");
	/* clobbers*/
	return ((unsigned long long) lo) | (((unsigned long long) hi) << 32);
}

ticks startRDTSC(void) {
	return rdtsc();
}

ticks stopRDTSCP(void) {
	return rdtsc();
}



#include <wmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

// Gueron fig. 5
__m128i gfmul_fig5(__m128i a, __m128i b) {
	__m128i tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
	tmp3 = _mm_clmulepi64_si128(a, b, 0x00);
	tmp4 = _mm_clmulepi64_si128(a, b, 0x10);
	tmp5 = _mm_clmulepi64_si128(a, b, 0x01);
	tmp6 = _mm_clmulepi64_si128(a, b, 0x11);

	tmp4 = _mm_xor_si128(tmp4, tmp5);
	tmp5 = _mm_slli_si128(tmp4, 8);
	tmp4 = _mm_srli_si128(tmp4, 8);
	tmp3 = _mm_xor_si128(tmp3, tmp5);
	tmp6 = _mm_xor_si128(tmp6, tmp4);

	tmp7 = _mm_srli_epi32(tmp3, 31);
	tmp8 = _mm_srli_epi32(tmp6, 31);
	tmp3 = _mm_slli_epi32(tmp3, 1);
	tmp6 = _mm_slli_epi32(tmp6, 1);

	tmp9 = _mm_srli_si128(tmp7, 12);
	tmp8 = _mm_slli_si128(tmp8, 4);
	tmp7 = _mm_slli_si128(tmp7, 4);
	tmp3 = _mm_or_si128(tmp3, tmp7);
	tmp6 = _mm_or_si128(tmp6, tmp8);
	tmp6 = _mm_or_si128(tmp6, tmp9);

	tmp7 = _mm_slli_epi32(tmp3, 31);
	tmp8 = _mm_slli_epi32(tmp3, 30);
	tmp9 = _mm_slli_epi32(tmp3, 25);

	tmp7 = _mm_xor_si128(tmp7, tmp8);
	tmp7 = _mm_xor_si128(tmp7, tmp9);
	tmp8 = _mm_srli_si128(tmp7, 4);
	tmp7 = _mm_slli_si128(tmp7, 12);
	tmp3 = _mm_xor_si128(tmp3, tmp7);

	tmp2 = _mm_srli_epi32(tmp3, 1);
	tmp4 = _mm_srli_epi32(tmp3, 2);
	tmp5 = _mm_srli_epi32(tmp3, 7);
	tmp2 = _mm_xor_si128(tmp2, tmp4);
	tmp2 = _mm_xor_si128(tmp2, tmp5);
	tmp2 = _mm_xor_si128(tmp2, tmp8);
	tmp3 = _mm_xor_si128(tmp3, tmp2);
	tmp6 = _mm_xor_si128(tmp6, tmp3);
	return tmp6;
}


// Gueron fig. 8
__m128i gfmul_fig8(__m128i H1, __m128i H2, __m128i H3, __m128i H4, __m128i X1, __m128i X2, __m128i X3, __m128i X4 )
{
	__m128i H1_X1_lo, H1_X1_hi,H2_X2_lo, H2_X2_hi, H3_X3_lo, H3_X3_hi, H4_X4_lo, H4_X4_hi, lo, hi;
    __m128i tmp0, tmp1, tmp2, tmp3; __m128i tmp4, tmp5, tmp6, tmp7; __m128i tmp8, tmp9;
    H1_X1_lo = _mm_clmulepi64_si128(H1, X1,0x00);
    H2_X2_lo = _mm_clmulepi64_si128(H2, X2,0x00);
    H3_X3_lo = _mm_clmulepi64_si128(H3, X3, 0x00);
    H4_X4_lo = _mm_clmulepi64_si128(H4, X4,0x00);
    lo = _mm_xor_si128(H1_X1_lo, H2_X2_lo);
    lo = _mm_xor_si128(lo, H3_X3_lo);
    lo = _mm_xor_si128(lo, H4_X4_lo);
    H1_X1_hi = _mm_clmulepi64_si128(H1, X1, 0x11);
    H2_X2_hi = _mm_clmulepi64_si128(H2, X2, 0x11);
    H3_X3_hi = _mm_clmulepi64_si128(H3, X3, 0x11);
    H4_X4_hi = _mm_clmulepi64_si128(H4, X4,0x11);
    hi = _mm_xor_si128(H1_X1_hi, H2_X2_hi);
    hi = _mm_xor_si128(hi, H3_X3_hi);
    hi = _mm_xor_si128(hi, H4_X4_hi);
    tmp0 = _mm_shuffle_epi32(H1, 78);
    tmp4 = _mm_shuffle_epi32(X1, 78);
    tmp0 = _mm_xor_si128(tmp0, H1);
    tmp4 = _mm_xor_si128(tmp4, X1);
    tmp1 = _mm_shuffle_epi32(H2, 78);
    tmp5 = _mm_shuffle_epi32(X2, 78);
    tmp1 = _mm_xor_si128(tmp1, H2);
    tmp5 = _mm_xor_si128(tmp5, X2);
    tmp2 = _mm_shuffle_epi32(H3, 78);

    tmp6 = _mm_shuffle_epi32(X3, 78);
    tmp2 = _mm_xor_si128(tmp2, H3);
    tmp6 = _mm_xor_si128(tmp6, X3);
    tmp3 = _mm_shuffle_epi32(H4, 78);
    tmp7 = _mm_shuffle_epi32(X4, 78);
    tmp3 = _mm_xor_si128(tmp3, H4);
    tmp7 = _mm_xor_si128(tmp7, X4);
    tmp0 = _mm_clmulepi64_si128(tmp0,tmp4, 0x00);
    tmp1 = _mm_clmulepi64_si128(tmp1,tmp5, 0x00);
    tmp2 = _mm_clmulepi64_si128(tmp2,tmp6, 0x00);
    tmp3 = _mm_clmulepi64_si128(tmp3,tmp7, 0x00);
    tmp0 = _mm_xor_si128(tmp0, lo);
    tmp0 = _mm_xor_si128(tmp0, hi);
    tmp0 = _mm_xor_si128(tmp1, tmp0);
    tmp0 = _mm_xor_si128(tmp2, tmp0);
    tmp0 = _mm_xor_si128(tmp3, tmp0);
    tmp4 = _mm_slli_si128(tmp0, 8);
    tmp0 = _mm_srli_si128(tmp0, 8);
    lo = _mm_xor_si128(tmp4, lo);
    hi = _mm_xor_si128(tmp0, hi);
    tmp3 = lo;
    tmp6 = hi;
    tmp7 = _mm_srli_epi32(tmp3, 31);
    tmp8 = _mm_srli_epi32(tmp6, 31);
    tmp3 = _mm_slli_epi32(tmp3, 1);
    tmp6 = _mm_slli_epi32(tmp6, 1);
    tmp9 = _mm_srli_si128(tmp7, 12);
    tmp8 = _mm_slli_si128(tmp8, 4);
    tmp7 = _mm_slli_si128(tmp7, 4);
    tmp3 = _mm_or_si128(tmp3, tmp7);
    tmp6 = _mm_or_si128(tmp6, tmp8);
    tmp6 = _mm_or_si128(tmp6, tmp9);
    tmp7 = _mm_slli_epi32(tmp3, 31);
    tmp8 = _mm_slli_epi32(tmp3, 30);
    tmp9 = _mm_slli_epi32(tmp3, 25);
    tmp7 = _mm_xor_si128(tmp7, tmp8);
    tmp7 = _mm_xor_si128(tmp7, tmp9);
    tmp8 = _mm_srli_si128(tmp7, 4);
    tmp7 = _mm_slli_si128(tmp7, 12);
    tmp3 = _mm_xor_si128(tmp3, tmp7);
    tmp2 = _mm_srli_epi32(tmp3, 1);
    tmp4 = _mm_srli_epi32(tmp3, 2);
    tmp5 = _mm_srli_epi32(tmp3, 7);
    tmp2 = _mm_xor_si128(tmp2, tmp4);
    tmp2 = _mm_xor_si128(tmp2, tmp5);
    tmp2 = _mm_xor_si128(tmp2, tmp8);
    tmp3 = _mm_xor_si128(tmp3, tmp2);
    tmp6 = _mm_xor_si128(tmp6, tmp3);
    return tmp6;
}

void printusage(char * command) {
	printf(" Usage: %s ", command);
}

void printme32(__m128i v1) {
	printf(" %u %u %u %u  ", _mm_extract_epi32(v1,0), _mm_extract_epi32(v1,1), _mm_extract_epi32(v1,2), _mm_extract_epi32(v1,3));
}

// inefficient but we do not care
int equal(__m128i a, __m128i b) {
	return (_mm_extract_epi32(a,0) == _mm_extract_epi32(b,0))
			&& (_mm_extract_epi32(a,1) == _mm_extract_epi32(b,1))
			&& (_mm_extract_epi32(a,2) == _mm_extract_epi32(b,2))
			&& (_mm_extract_epi32(a,3) == _mm_extract_epi32(b,3));
}

int main(int argc, char ** arg) {
	int SHORTTRIALS = 1000000;// should be divisible by 4
	assert(SHORTTRIALS/4*4 == SHORTTRIALS);
	int HowManyRepeats = 5;
	int elapsed1;
	int j,k;

	ticks bef, aft;
	struct timeval start, finish;
	int c;
	while ((c = getopt(argc, arg, "h")) != -1)
		switch (c) {
		case 'h':
			printusage(arg[0]);
			return 0;
		default:
			abort();
		}
	__m128i A;
	__m128i B;
	//////////////////////
	// A is our polynomial coefficient, and the string is made of Bs
	// so we compute A^n * B + ... A^2 * B + A*B.
	// In real-life, you would need to load different values for B but
	// this can only add some overhead that reduces the difference between
	// schemes.
	//////////////////////

	for (k = 0; k < HowManyRepeats; ++k) {
		A = _mm_set1_epi32(~k);
		B = _mm_set1_epi32(~(k+1));

		printf("test #%d  \n", k + 1);
		//
		gettimeofday(&start, 0);
		bef = startRDTSC();
		__m128i X1 = _mm_setzero_si128();
		for (j = 0; j < SHORTTRIALS; ++j) {
			X1 = _mm_xor_si128(X1,B);
			X1 = mul128by128to128_lazymod127(A, X1);
		}
		aft = stopRDTSCP();
		gettimeofday(&finish, 0);
		elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
				+ (finish.tv_usec - start.tv_usec));

		printf(
				"[basic 127-bit technique] CPU cycle/mult = %f \t billions of mult per second =  %f    \n",
				(aft - bef) * 1.0 / (4.0 * SHORTTRIALS ),
				(4.0 * SHORTTRIALS ) / (1000. * elapsed1));
		//
		gettimeofday(&start, 0);
		bef = startRDTSC();
		__m128i X2 = _mm_setzero_si128();
		for (j = 0; j < SHORTTRIALS; ++j) {
			X2 = _mm_xor_si128(X2,B);
			X2 = gfmul_fig5(A, X2);
		}
		aft = stopRDTSCP();
		gettimeofday(&finish, 0);
		elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
				+ (finish.tv_usec - start.tv_usec));

		printf(
				"[basic Gueron GHASH technique] CPU cycle/mult = %f \t billions of mult per second =  %f    \n",
				(aft - bef) * 1.0 / (4.0 * SHORTTRIALS ),
				(4.0 * SHORTTRIALS ) / (1000. * elapsed1));
		//
		__m128i A2 = gfmul_fig5(A, A);
		__m128i A3 = gfmul_fig5(A2, A);
		__m128i A4 = gfmul_fig5(A2, A2);
		gettimeofday(&start, 0);
		bef = startRDTSC();
		__m128i X3 = _mm_setzero_si128();;
		for (j = 0; j < SHORTTRIALS/4; ++j) {
			X3 = _mm_xor_si128(X3,B);
			X3 = gfmul_fig8(A, A2, A3, A4,  B, B, B, X3);
		}
		aft = stopRDTSCP();
		gettimeofday(&finish, 0);
		elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
				+ (finish.tv_usec - start.tv_usec));

		printf(
				"[aggregated(4) Gueron GHASH technique] CPU cycle/mult = %f \t billions of mult per second =  %f    \n",
				(aft - bef) * 1.0 / (4.0 * SHORTTRIALS ),
				(4.0 * SHORTTRIALS ) / (1000. * elapsed1));

		//
		printf("X1=");
		printme32(X1);
		printf("X2=");
		printme32(X2);
		printf("X3=");
		printme32(X3);
		assert(equal(X2,X3));
		printf("\n");
		printf("\n");
	}

	printf("\n");
}




