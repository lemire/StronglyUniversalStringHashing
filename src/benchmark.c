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
// Gabriele Paoloni, How to Benchmark Code Execution Times on Intel® IA-32 and IA-64 Instruction Set Architectures
// September 2010
// http://edc.intel.com/Link.aspx?id=3954

static __inline__ ticks fancystartRDTSC(void) {
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
}

#include "hashfunctions32bits.h"
#include "hashfunctions64bits.h"

#ifdef __PCLMUL__

#include "clmulhashfunctions32bits.h"
#include "clmulhashfunctions64bits.h"
#include "clmulpoly64bits.h"

#define HowManyFunctions 12
#define HowManyFunctions64 16


hashFunction64 funcArr64[HowManyFunctions64] = {&hashCity,
		&hashMMH_NonPyramidal,
		&hashNH64,
		&hashGaloisFieldfast64_precomp,
	&hashGaloisFieldfast64halfunrolled,
	&hashGaloisFieldfast64halfunrolled_precomp,
	&hashGaloisFieldPoly64,
	&precomphashGaloisFieldPoly64,&fasthashGaloisFieldPoly64_2_noprecomp,
	&fasthashGaloisFieldPoly64_2,&fasthashGaloisFieldPoly64_4,&fasthashGaloisFieldPoly64_8,
	&fasthashGaloisFieldPoly64_16,&halfhashGaloisFieldPoly64_8,&halfhashGaloisFieldPoly64_16,
	&clmulgarbage
};

hashFunction funcArr[HowManyFunctions] = {&hashGaloisFieldMultilinear,
	&hashGaloisFieldMultilinearHalfMultiplications, &hashMultilinear,&hashMultilinear2by2 ,
	&hashMultilinearhalf, &hashMultilineardouble,
	&hashNH,&hashRabinKarp, &hashFNV1, &hashFNV1a, &hashSAX,&pyramidal_Multilinear
};

const char* functionnames64[HowManyFunctions64] = {
	"Google's City                       ",
	"Non-pyramidal MMH                   ",
	"Simple NH64 (like VMAC)             ",
	"GFMultilinear (64-bit regular pre)  ",
	"GFMultilinear (64-bit half, unrol)  ",
	"GFMultilinear(64-bit half,unrol,pre)",
	"hashGaloisFieldPoly64               ",
	"hashGaloisFieldPoly64  (precomp)    ",
	"fasthashGaloisFieldPoly64 (nopre 2) ",
	"fasthashGaloisFieldPoly64 (2)       ",
	"fasthashGaloisFieldPoly64 (4)       ",
	"fasthashGaloisFieldPoly64 (8)       ",
	"fasthashGaloisFieldPoly64 (16)      ",
	"halfhashGaloisFieldPoly64 (8)       ",
	"halfhashGaloisFieldPoly64 (16)      ",
	"garbage                             ",

};

const char* functionnames[HowManyFunctions] = {
	"GFMultilinear   (strongly universal)",
	"GFMultilinearhalf   (str. universal)",
	"Multilinear     (strongly universal)",
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
#else

#define HowManyFunctions 10
#define HowManyFunctions64 3

hashFunction funcArr[HowManyFunctions] = { &hashMultilinear,
		&hashMultilinear2by2, &hashMultilinearhalf, &hashMultilineardouble,
		&hashNH, &hashRabinKarp, &hashFNV1, &hashFNV1a, &hashSAX,
		&pyramidal_Multilinear };
hashFunction64 funcArr64[HowManyFunctions64] = {&hashCity, &hashMMH_NonPyramidal, &hashNH64 };

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
		"Pyramidal multilinear (a. univ.)    " };
const char* functionnames64[HowManyFunctions64] = {
		"Google's City                       ",
		"Non-pyramidal MMH                   ",
		"Simple NH64 (like VMAC)             ", };

#endif

int main(int c, char ** arg) {
	(void) (c);
	(void) (arg);
	const int N = 1024; // should be divisible by two!
	const int SHORTTRIALS = 1000000;
	const int HowManyRepeats = 3;
	int i, k, j;
	int elapsed;
	hashFunction thisfunc;
	const char * functionname;
	ticks bef, aft;
	struct timeval start, finish;
	uint64_t randbuffer[N + 3] __attribute__ ((aligned (16)));
	uint32_t sumToFoolCompiler = 0;
	uint32_t intstring[N] __attribute__ ((aligned (16))); // // could force 16-byte alignment with  __attribute__ ((aligned (16)));
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
		printf("test #%d (64-bit hash values) ", k + 1);
		printf("(%d bytes) \n", N * 4);

		hashFunction64 thisfunc64;
		for (i = 0; i < HowManyFunctions64; ++i) {
			sumToFoolCompiler = 0;
			thisfunc64 = funcArr64[i];
			functionname = functionnames64[i];
			printf("%s ", functionname);
			fflush(stdout);
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
			printf("# ignore this #%d\n", sumToFoolCompiler);

		}
		printf("\n");
		printf("test #%d (32-bit hash values)\n", k + 1);
		for (i = 0; i < HowManyFunctions; ++i) {
			sumToFoolCompiler = 0;
			thisfunc = funcArr[i];
			functionname = functionnames[i];
			printf("%s ", functionname);
			fflush(stdout);
			gettimeofday(&start, 0);
			bef = startRDTSC();
			for (j = 0; j < SHORTTRIALS; ++j)
				sumToFoolCompiler += thisfunc(&randbuffer[0], &intstring[0], N);
			aft = stopRDTSCP();
			gettimeofday(&finish, 0);
			elapsed = (1000000 * (finish.tv_sec - start.tv_sec)
					+ (finish.tv_usec - start.tv_usec));
			printf(
					"CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",
					(aft - bef) * 1.0 / (4.0 * SHORTTRIALS * N),
					(4.0 * SHORTTRIALS * N) / (1000. * elapsed));
			printf("# ignore this #%d\n", sumToFoolCompiler);

		}
		printf("\n");
	}
	printf("# ignore this #%d\n", sumToFoolCompiler);

}

