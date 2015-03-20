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
	int SHORTTRIALS = 100000;// should be divisible by 4
	assert(SHORTTRIALS/4*4 == SHORTTRIALS);
	int N = SHORTTRIALS * sizeof(__m128i)/sizeof(uint32_t);
	uint32_t intstring[N] __attribute__ ((aligned (16))); // // could force 16-byte alignment with  __attribute__ ((aligned (16)));
	int HowManyRepeats = 5;
	int elapsed1;
	int i, j,k;

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
	for (i = 0; i < N; ++i) {
		intstring[i] = rand();
	}
	__m128i * input = (__m128i *) intstring;
	__m128i A;
	//////////////////////
	// A is our polynomial coefficient, and the string is made of Bs
	// so we compute A^n * Bn + ... A^2 * B2 + A*B1.
	//////////////////////

	for (k = 0; k < HowManyRepeats; ++k) {
		A = _mm_set1_epi32(~k);
		// next we see the two highest bits to zero to get a 126-bit key
		A = _mm_and_si128(A,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x3fffffff));// setting two highest bits to zero

		printf("test #%d using %lu kB of random input data \n", k + 1 , N*sizeof(uint32_t)/1024);
		//
		gettimeofday(&start, 0);
		bef = startRDTSC();
		__m128i X1 = _mm_setzero_si128();
		for (j = 0; j < SHORTTRIALS; ++j) {
			__m128i B = _mm_load_si128(input + j );
			X1 = _mm_xor_si128(X1,B);
			X1 = mul128by128to128_lazymod127(A, X1);
		}
		aft = stopRDTSCP();
		gettimeofday(&finish, 0);
		elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
				+ (finish.tv_usec - start.tv_usec));

		printf(
				"[basic 126-bit technique] CPU cycle/mult = %f \t billions of mult per second =  %f    \n",
				(aft - bef) * 1.0 / (4.0 * SHORTTRIALS ),
				(4.0 * SHORTTRIALS ) / (1000. * elapsed1));
		//
		gettimeofday(&start, 0);
		bef = startRDTSC();
		__m128i X2 = _mm_setzero_si128();
		for (j = 0; j < SHORTTRIALS; ++j) {
			__m128i B = _mm_load_si128(input + j );
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
			__m128i B1 = _mm_load_si128(input + 4 * j );
			__m128i B2 = _mm_load_si128(input + 4 * j + 1);
			__m128i B3 = _mm_load_si128(input + 4 * j + 2);
			__m128i B4 = _mm_load_si128(input + 4 * j + 3);
			X3 = _mm_xor_si128(X3,B1);
			X3 = gfmul_fig8(A, A2, A3, A4,  B4, B3, B2, X3);
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
		__m128i A2lazy = mul128by128to128_lazymod127(A, A);
		gettimeofday(&start, 0);
		bef = startRDTSC();
		__m128i X4 = _mm_setzero_si128();;
		for (j = 0; j < SHORTTRIALS/2; ++j) {
			__m128i B1 = _mm_load_si128(input + 2 * j );
			__m128i B2 = _mm_load_si128(input + 2 * j + 1);
			X4 = _mm_xor_si128(X4,B1);
			X4 = mul128by128to128_lazymod127_2by2(A, A2lazy, B2, X4);
		}
		aft = stopRDTSCP();
		gettimeofday(&finish, 0);
		elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
					+ (finish.tv_usec - start.tv_usec));

		printf(
					"[aggregated(2) 126-bit technique] CPU cycle/mult = %f \t billions of mult per second =  %f    \n",
					(aft - bef) * 1.0 / (4.0 * SHORTTRIALS ),
					(4.0 * SHORTTRIALS ) / (1000. * elapsed1));

		//
		//
		__m128i A3lazy = mul128by128to128_lazymod127(A2lazy, A);
		__m128i A4lazy = mul128by128to128_lazymod127(A2lazy, A2lazy);

		gettimeofday(&start, 0);
		bef = startRDTSC();
		__m128i X5 = _mm_setzero_si128();
		for (j = 0; j < SHORTTRIALS/4; ++j) {
			__m128i B1 = _mm_load_si128(input + 4 * j );
			__m128i B2 = _mm_load_si128(input + 4 * j + 1);
			__m128i B3 = _mm_load_si128(input + 4 * j + 2);
			__m128i B4 = _mm_load_si128(input + 4 * j + 3);
			X5 = _mm_xor_si128(X5,B1);
			X5 = mul128by128to128_lazymod127_4by4(A, A2lazy, A3lazy, A4lazy, B4, B3, B2, X5);
		}
		aft = stopRDTSCP();
		gettimeofday(&finish, 0);
		elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
					+ (finish.tv_usec - start.tv_usec));

		printf(
					"[aggregated(4) 126-bit technique] CPU cycle/mult = %f \t billions of mult per second =  %f    \n",
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
		printf("X4=");
		printme32(X4);
		assert(equal(X1,X4));// they could differ in theory...
		printf("X5=");
		printme32(X5);
		assert(equal(X1,X5));// they could differ in theory...
		printf("\n");
		printf("\n");
	}

	printf("\n");
}




