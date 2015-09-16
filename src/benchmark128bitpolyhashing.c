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
#include "ghash.h"


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





void printusage(char * command) {
    printf(" Usage: %s ", command);
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




