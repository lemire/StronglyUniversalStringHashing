/*
 * Benchmarking carry-less multiplications 128x128 to 256
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


void mul128by128to256( __m128i A, __m128i B, __m128i * Alow,  __m128i * Ahigh) {
    __m128i Alowtmp, Ahightmp, Amix1, Amix2, Amix;
    Amix1 = _mm_clmulepi64_si128(A,B,0x01);
    Amix2 = _mm_clmulepi64_si128(A,B,0x10);
    Alowtmp = _mm_clmulepi64_si128(A,B,0x00);
    Ahightmp = _mm_clmulepi64_si128(A,B,0x11);
    Amix = _mm_xor_si128(Amix1,Amix2);
    Amix1 = _mm_slli_si128(Amix,8);
    Amix2 = _mm_srli_si128(Amix,8);
    *Alow = _mm_xor_si128(Alowtmp,Amix1);
    *Ahigh = _mm_xor_si128(Ahightmp,Amix2);
}

void mul128by128to256_gueron_fig5( __m128i a, __m128i b, __m128i * Alow,  __m128i * Ahigh) {
    __m128i tmp3, tmp4, tmp5, tmp6;
    tmp3 = _mm_clmulepi64_si128(a, b, 0x00);
    tmp4 = _mm_clmulepi64_si128(a, b, 0x10);
    tmp5 = _mm_clmulepi64_si128(a, b, 0x01);
    tmp6 = _mm_clmulepi64_si128(a, b, 0x11);
    tmp4 = _mm_xor_si128(tmp4, tmp5);
    tmp5 = _mm_slli_si128(tmp4, 8);
    tmp4 = _mm_srli_si128(tmp4, 8);
    *Alow = _mm_xor_si128(tmp3, tmp5);
    *Ahigh = _mm_xor_si128(tmp6, tmp4);
}

void mul128by128to256_gueron_fig7( __m128i a, __m128i b, __m128i * Alow,  __m128i * Ahigh) {
    __m128i tmp3, tmp4, tmp5, tmp6;
    tmp3 = _mm_clmulepi64_si128(a, b, 0x00);
    tmp6 = _mm_clmulepi64_si128(a, b, 0x11);
    tmp4 = _mm_shuffle_epi32(a,78);// 78 = 0b 01 00 11 10 ; 2->0, 3->1 0->2 1->3 // flips high and low 64-bit words
    tmp5 = _mm_shuffle_epi32(b,78);
    tmp4 = _mm_xor_si128(tmp4, a);
    tmp5 = _mm_xor_si128(tmp5, b);
    tmp4 = _mm_clmulepi64_si128(tmp4, tmp5,0x00);// (A1+A0)*(B1+B0) = A1*B1 + A1*B0 + A0*B0 + A0*B0
    tmp4 = _mm_xor_si128(tmp4, tmp3); //(A1*B1 + A1*B0 + A0*B1 + A0*B0) - A0*B0
    tmp4 = _mm_xor_si128(tmp4, tmp6); //(A1*B1 + A1*B0 + A0*B1 ) - A1*B1 = ( A1*B0 + A0*B1 )
    tmp5 = _mm_slli_si128(tmp4, 8);
    tmp4 = _mm_srli_si128(tmp4, 8);
    *Alow = _mm_xor_si128(tmp3, tmp5);
    *Ahigh = _mm_xor_si128(tmp6, tmp4);
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
    int SHORTTRIALS = 1000000;
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
    for (k = 0; k < HowManyRepeats; ++k) {
        A = _mm_set1_epi32(~k);
        B = _mm_set1_epi32(~(k+1));
        printf("test #%d \n ", k + 1);
        //
        gettimeofday(&start, 0);
        bef = startRDTSC();
        __m128i Alow1 = A;
        __m128i Ahigh1 = B;
        for (j = 0; j < SHORTTRIALS; ++j) {
            mul128by128to256(Alow1,Ahigh1,&Alow1,&Ahigh1);
            Alow1 = _mm_xor_si128(Alow1,A);
            Ahigh1 = _mm_xor_si128(Ahigh1,B);
        }
        aft = stopRDTSCP();
        gettimeofday(&finish, 0);
        elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
                    + (finish.tv_usec - start.tv_usec));

        printf(
            "[basic technique] CPU cycle/mult = %f \t billions of mult per second =  %f    \n",
            (aft - bef) * 1.0 / (4.0 * SHORTTRIALS ),
            (4.0 * SHORTTRIALS ) / (1000. * elapsed1));
        //
        gettimeofday(&start, 0);
        bef = startRDTSC();
        __m128i Alow2 = A;
        __m128i Ahigh2 = B;
        for (j = 0; j < SHORTTRIALS; ++j) {
            mul128by128to256_gueron_fig5(Alow2,Ahigh2,&Alow2,&Ahigh2);
            Alow2 = _mm_xor_si128(Alow2,A);
            Ahigh2 = _mm_xor_si128(Ahigh2,B);
        }
        aft = stopRDTSCP();
        gettimeofday(&finish, 0);
        elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
                    + (finish.tv_usec - start.tv_usec));

        printf(
            "[basic Gueron technique] CPU cycle/mult = %f \t billions of mult per second =  %f    \n",
            (aft - bef) * 1.0 / (4.0 * SHORTTRIALS ),
            (4.0 * SHORTTRIALS ) / (1000. * elapsed1));
        //
        gettimeofday(&start, 0);
        bef = startRDTSC();
        __m128i Alow3 = A;
        __m128i Ahigh3 = B;
        for (j = 0; j < SHORTTRIALS; ++j) {
            mul128by128to256_gueron_fig7(Alow3,Ahigh3,&Alow3,&Ahigh3);
            Alow3 = _mm_xor_si128(Alow3,A);
            Ahigh3 = _mm_xor_si128(Ahigh3,B);
        }
        aft = stopRDTSCP();
        gettimeofday(&finish, 0);
        elapsed1 = (1000000 * (finish.tv_sec - start.tv_sec)
                    + (finish.tv_usec - start.tv_usec));

        printf(
            "[3-mult Gueron technique] CPU cycle/mult = %f \t billions of mult per second =  %f    \n",
            (aft - bef) * 1.0 / (4.0 * SHORTTRIALS ),
            (4.0 * SHORTTRIALS ) / (1000. * elapsed1));


        printme32(Alow1);
        printme32(Alow2);
        printme32(Alow3);
        assert( equal(Alow1,Alow2) );
        assert( equal(Alow2,Alow3) );
        printf("\n");
        printme32(Ahigh1);
        printme32(Ahigh2);
        printme32(Ahigh3);
        assert( equal(Ahigh1,Ahigh2) );
        assert( equal(Ahigh1,Ahigh3) );
        printf("\n");
        printf("\n");
    }

    printf("\n");
}

