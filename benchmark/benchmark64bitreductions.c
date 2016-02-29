
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>
#ifdef __AVX__
#define __PCLMUL__ 1
#endif

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

#include "clmul.h"

void force_computation(uint32_t forcedValue) {
    // make sure forcedValue has to be computed, but avoid output (unless unlucky)
    if (forcedValue % 277387 == 17)
        printf("wow, what a coincidence! (in benchmark.c)");
}

void printusage(char * command) {
    printf(" Usage: %s ", command);
}

// WARNING: HIGH 64 BITS CONTAIN GARBAGE, must call _mm_cvtsi128_si64 to get
// meaningful bits.
__m128i barrettWithoutPrecomputation64_si128(__m128i A) {
    ///http://www.jjj.de/mathdata/minweight-primpoly.txt
    // it is important, for the algo. we have chosen that 4 is smaller
    // equal than 32=64/2

    //const int n = 64;	// degree of the polynomial
    //const __m128i C = _mm_set_epi64x(1U,
    //		(1U << 4) + (1U << 3) + (1U << 1) + (1U << 0));	// C is the irreducible poly. (64,4,3,1,0)
    const __m128i C = _mm_cvtsi64_si128((1U << 4) + (1U << 3) + (1U << 1) + (1U << 0));
    /////////////////
    /// This algo. requires two multiplications (_mm_clmulepi64_si128)
    /// They are probably the bottleneck.
    /// Note: Barrett's original algorithm also required two multiplications.
    ////////////////
    //assert(n / 8 == 8);
    __m128i Q2 = _mm_clmulepi64_si128(A, C, 0x01);
    Q2 = _mm_xor_si128(Q2, A);
    const __m128i Q4 = _mm_clmulepi64_si128(Q2, C, 0x01);
    const __m128i final = _mm_xor_si128(A, Q4);
    return final;	/// WARNING: HIGH 64 BITS CONTAIN GARBAGE
}

uint64_t barrettWithoutPrecomputation64(__m128i A) {
    const __m128i final = barrettWithoutPrecomputation64_si128(A);
    return _mm_cvtsi128_si64(final);
}

int main(int argc, char ** arg) {
    int N = 1024;
    int SHORTTRIALS = 100000;
    int HowManyRepeats = 3;
    int elapsed1, elapsed2;
    int i,j,k;
    int sumToFoolCompiler1, sumToFoolCompiler2;

    ticks bef, aft;
    struct timeval start, finish;
    uint32_t intstring[N] __attribute__ ((aligned (16))); // // could force 16-byte alignment with  __attribute__ ((aligned (16)));
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

    printf(
        "Reporting the number of cycles per byte and the billions of bytes processed per second.\n");
    for (k = 0; k < HowManyRepeats; ++k) {
        printf("test #%d (reduction) ", k + 1);
        printf("(%d bytes) \n", N * 4);
        sumToFoolCompiler1 = 0;
        sumToFoolCompiler2 = 0;
        __m128i * data = (__m128i *) &intstring[0];
        gettimeofday(&start, 0);
        bef = startRDTSC();
        for (j = 0; j < SHORTTRIALS; ++j)
            for (i = 0; i < N / 4; ++i)
                sumToFoolCompiler1 += precompReduction64(data[i]);
        aft = stopRDTSCP();
        gettimeofday(&finish, 0);
        elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
                    + (finish.tv_usec - start.tv_usec));

        printf(
            "[fast technique] CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",
            (aft - bef) * 1.0 / (4.0 * SHORTTRIALS * N),
            (4.0 * SHORTTRIALS * N) / (1000. * elapsed1));
        force_computation (sumToFoolCompiler1);
        gettimeofday(&start, 0);
        bef = startRDTSC();
        for (j = 0; j < SHORTTRIALS; ++j)
            for (i = 0; i < N / 4; ++i)
                sumToFoolCompiler2 += barrettWithoutPrecomputation64(data[i]);
        aft = stopRDTSCP();
        gettimeofday(&finish, 0);
        elapsed2 = (1000000 * (finish.tv_sec - start.tv_sec)
                    + (finish.tv_usec - start.tv_usec));
        printf(
            "[noprecomp     ] CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",
            (aft - bef) * 1.0 / (4.0 * SHORTTRIALS * N),
            (4.0 * SHORTTRIALS * N) / (1000. * elapsed2));
        printf("speed ratio = %f \n",1.*elapsed2/elapsed1);
        force_computation (sumToFoolCompiler2);
    }
    printf("\n");
}

