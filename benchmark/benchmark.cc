/////////////////////////////////////
// This C code is a companion to the paper
//
// Reference: Owen Kaser and Daniel Lemire, Strongly universal string hashing is fast, Computer Journal
// http://arxiv.org/abs/1202.4961
//
// It shows that we can compute strongly universal hash functions very quickly.
/////////////////////////////////////

//
// this code will hash strings of 32-bit characters. To use on
// strings of 8-bit characters, you may need some adequate padding.
//
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
// start and stop are as recommended by
// Gabriele Paoloni, How to Benchmark Code Execution Times on Intel IA-32 and IA-64 Instruction Set Architectures
// September 2010
// http://edc.intel.com/Link.aspx?id=3954

/*static __inline__ ticks fancystartRDTSC(void) {
	unsigned cycles_low, cycles_high;
	asm volatile ("CPUID\n\t"
			"RDTSC\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t": "=r" (cycles_high), "=r" (cycles_low)::
			"%rax", "%rbx", "%rcx", "%rdx");
	return ((ticks) cycles_high << 32) | cycles_low;
}

static __inline__ ticks fancystopRDTSCP(void) {
	unsigned cycles_low, cycles_high;
/// This should work fine on most machines, if the RDTSCP thing
/// fails for you, use the  rdtsc() call instead.
	asm volatile("RDTSCP\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t"
			"CPUID\n\t": "=r" (cycles_high), "=r" (cycles_low):: "%rax",
			"%rbx", "%rcx", "%rdx");
	return ((ticks) cycles_high << 32) | cycles_low;
}*/

extern "C" {

#include "hashfunctions32bits.h"
#include "hashfunctions64bits.h"

}

#ifdef __PCLMUL__

extern "C" {
#include "clmulhashfunctions32bits.h"
#include "clmulhashfunctions64bits.h"
#include "clmulpoly64bits.h"
#include "ghash.h"
#include "clmulhierarchical64bits.h"
}

#include "treehash/binary-treehash.hh"
#include "treehash/generic-treehash.hh"


#define HowManyFunctions 13
#define HowManyFunctions64 11

hashFunction64 funcArr64[HowManyFunctions64] = {&hashCity,
                                                &hashVHASH64,
                                                &CLHASH,
                                                &hashGaloisFieldfast64_precomp_unroll,
                                                &hashGaloisFieldfast64halfunrolled_precomp,
                                                &hashSipHash,
                                                &GHASH64bit,
                                                &generic_treehash<BoostedZeroCopyGenericBinaryTreehash, CLNH, 7>,
                                                &generic_treehash<BoostedZeroCopyGenericBinaryTreehash, NHCL, 7>,
                                                &generic_treehash<BoostedZeroCopyGenericBinaryTreehash, NH, 7>,
                                                &hashPMP64,
                                               };

hashFunction funcArr[HowManyFunctions] = {&hashGaloisFieldMultilinear,
    &hashGaloisFieldMultilinearHalfMultiplications, &hashMultilinear,
    &hashMultilinear2by2, &hashMultilinearhalf, &hashMultilineardouble, &hashNH,
    &hashRabinKarp, &hashFNV1, &hashFNV1a, &hashSAX, &pyramidal_Multilinear, &pdp32avx};

const char* functionnames64[HowManyFunctions64] = {
    "Google's City                       ",
    "64-bit VHASH                        ",
    "64-bit CLHASH                       ",
    "GFMultilinear                       ",
    "GFMultilinear (half multiplication) ",
    "SipHash                             ",
    "GHASH                               ",
    "generic_tree<Boosted..., CLNH, 7>   ",
    "generic_tree<Boosted..., NHCL, 7>   ",
    "generic_tree<Boosted..., NH, 7>     ",
    "PMP64                               "
};

const char* functionnames[HowManyFunctions] = {
    "GFMultilinear   (poorly optimized)  ",
    "GFMultilinearhalf   (optimized)     ",
    "Multilinear     (strongly universal)",
    "Multilinear2x2  (strongly universal)",
    "Multilinearhalf (strongly universal)",
    "Multilineardouble (strongly u.)     ",
    "NH (high bits)                      ",
    "RabinKarp                           ",
    "FNV1                                ",
    "FNV1a                               ",
    "SAX                                 ",
    "Pyramidal multilinear (a. univ.)    ",
    "pdp32avx                            ",
};
#else

#define HowManyFunctions 10
#define HowManyFunctions64 3

hashFunction funcArr[HowManyFunctions] = { &hashMultilinear,
                                           &hashMultilinear2by2, &hashMultilinearhalf, &hashMultilineardouble,
                                           &hashNH, &hashRabinKarp, &hashFNV1, &hashFNV1a, &hashSAX,
                                           &pyramidal_Multilinear
                                         };
hashFunction64 funcArr64[HowManyFunctions64] = { &hashCity,
                                                 &hashMMH_NonPyramidal, &hashNH64
                                               };

const char* functionnames[HowManyFunctions] = {
    "Multilinear  (strongly universal)",
    "Multilinear2x2  (strongly universal)",
    "Multilinearhalf (strongly universal)",
    "Multilineardouble (strongly u.)     ",
    "NH (high bits)                      ",
    "RabinKarp                           ",
    "FNV1                                ",
    "FNV1a                               ",
    "SAX                                 ",
    "Pyramidal multilinear (a. univ.)    "
};
const char* functionnames64[HowManyFunctions64] = {
    "Google's City                       ",
    "Non-pyramidal MMH                   ",
    "Simple NH64 (like VMAC)             ",
};

#endif

void force_computation(uint32_t forcedValue) {
    // make sure forcedValue has to be computed, but avoid output (unless unlucky)
    if (forcedValue % 277387 == 17)
        printf("wow, what a coincidence! (in benchmark.c)");
    //printf("# ignore this #%d\n", force_computation);
}

void printusage(char * command) {
    printf(" Usage: %s -b (32|64)", command);
}

int main(int argc, char ** arg) {
    int N = 1024; // should be divisible by two!
    int SHORTTRIALS = 100000;
    int HowManyRepeats = 3;
    int bit = 64;
    int i, k, j;
    int elapsed;
    hashFunction thisfunc;
    const char * functionname;
    ticks bef, aft;
    struct timeval start, finish;
    uint64_t randbuffer[N + 3] __attribute__ ((aligned (32)));
    uint32_t sumToFoolCompiler = 0;
    uint32_t intstring[N] __attribute__ ((aligned (32)));
    int c;
    while ((c = getopt(argc, arg, "hb:")) != -1)
        switch (c) {
        case 'h':
            printusage(arg[0]);
            return 0;
        case 'b':
            bit = atoi(optarg);
            if ((bit != 32) && (bit != 64)) {
                printusage(arg[0]);
                return -1;
            }
            break;
        default:
            abort();
        }

    for (i = 0; i < N + 3; ++i) {
        randbuffer[i] = rand() | ((uint64_t)(rand()) << 32);
    }
    for (i = 0; i < N; ++i) {
        intstring[i] = rand();
    }
    printf(
        "For documentation, see Strongly universal string hashing is fast at http://arxiv.org/abs/1202.4961 \n");

    printf(
        "Reporting the number of cycles per byte and the billions of bytes processed per second.\n");
    for (k = 0; k < HowManyRepeats; ++k) {
        if (bit == 64) {
            printf("test #%d (64-bit hash values) ", k + 1);
            printf("(%d bytes) \n", N * 4);

            hashFunction64 thisfunc64;
            for (i = 0; i < HowManyFunctions64; ++i) {
                sumToFoolCompiler = 0;
                thisfunc64 = funcArr64[i];
                functionname = functionnames64[i];
                printf("%s ", functionname);
                fflush(stdout);
                sumToFoolCompiler += thisfunc64(&randbuffer[0],
                                                (uint64_t *) &intstring[0], N / 2);// we do not count the first run
                gettimeofday(&start, 0);
                bef = startRDTSC();
                assert(N / 2 * 2 == N);
                for (j = 0; j < SHORTTRIALS; ++j)
                    sumToFoolCompiler += thisfunc64(&randbuffer[0],
                                                    (uint64_t *) &intstring[0], N / 2);
                aft = stopRDTSCP();
                gettimeofday(&finish, 0);
                elapsed = (1000000 * (finish.tv_sec - start.tv_sec)
                           + (finish.tv_usec - start.tv_usec));
                printf(
                    "CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",
                    (aft - bef) * 1.0 / (4.0 * SHORTTRIALS * N),
                    (4.0 * SHORTTRIALS * N) / (1000. * elapsed));
                force_computation(sumToFoolCompiler);

            }
        } else {
            printf("test #%d (32-bit hash values)\n", k + 1);
            for (i = 0; i < HowManyFunctions; ++i) {
                sumToFoolCompiler = 0;
                thisfunc = funcArr[i];
                functionname = functionnames[i];
                printf("%s ", functionname);
                fflush(stdout);
                sumToFoolCompiler += thisfunc(&randbuffer[0], &intstring[0],
                                              N);// we do not count the first pass
                gettimeofday(&start, 0);
                bef = startRDTSC();
                for (j = 0; j < SHORTTRIALS; ++j)
                    sumToFoolCompiler += thisfunc(&randbuffer[0], &intstring[0],
                                                  N);
                aft = stopRDTSCP();
                gettimeofday(&finish, 0);
                elapsed = (1000000 * (finish.tv_sec - start.tv_sec)
                           + (finish.tv_usec - start.tv_usec));
                printf(
                    "CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",
                    (aft - bef) * 1.0 / (4.0 * SHORTTRIALS * N),
                    (4.0 * SHORTTRIALS * N) / (1000. * elapsed));
                force_computation(sumToFoolCompiler);

            }
        }
        printf("\n");
    }
    force_computation(sumToFoolCompiler); // printf("# ignore this #%d\n", sumToFoolCompiler);

}
