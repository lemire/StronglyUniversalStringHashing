
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
 #include <smmintrin.h>

//
// this is strongly universal. Condition: randomsource must be at least as long as 
// the string length  + 3.
//  
// Reference: Owen Kaser and Daniel Lemire, Strongly universal string hashing is fast, Computer Journal 
// http://arxiv.org/abs/1202.4961
//__attribute__ ((__target__ ("no-sse2"))) // GCC has buggy SSE2 code generation in some cases
uint32_t hashMultilinear(const uint64_t *  randomsource, const uint32_t *  string, const size_t length) {
    const uint32_t * const endstring = string + length;
    uint64_t sum = *(randomsource++);
    for(; string!= endstring; ++randomsource,++string ) {
        sum+= (*randomsource *  (uint64_t)(*string)) ;
    }
    sum += *randomsource;
    return (int) (sum>>32);
}



uint32_t hashMultilinearSSE41(const uint64_t *  randomsource, const uint32_t *  string, const size_t length) {
    uint64_t sum = *(randomsource++);
    randomsource++;// skipp
    __m128i counter =  _mm_setzero_si128 ();
    for(size_t i = 0; i < length;i+=4) {
    	__m128i randvec1 = _mm_load_si128 ((__m128i*) randomsource);
    	randomsource += 2;
    	__m128i randvec2 = _mm_load_si128 ((__m128i*) randomsource);
    	randomsource += 2;
    	__m128i datavec = _mm_load_si128 ((__m128i*) string);
    	string += 4;
    	
    	__m128i iterm1 =  _mm_slli_epi64 ( _mm_mul_epu32 (randvec1, datavec), 32);
    	counter = _mm_add_epi64 (counter, iterm1);
    	__m128i iterm2 =    _mm_mul_epu32 ( _mm_slli_si128(randvec1,4), datavec);
    	counter = _mm_add_epi64 (counter, iterm2);
    	datavec = _mm_slli_si128(datavec,4);
    	__m128i iterm3 =  _mm_slli_epi64 ( _mm_mul_epu32 (randvec2, datavec), 32);
    	counter = _mm_add_epi64 (counter, iterm3);
    	__m128i iterm4 =    _mm_mul_epu32 ( _mm_slli_si128(randvec2,4), datavec);
    	counter = _mm_add_epi64 (counter, iterm4);
    }

    sum +=  _mm_extract_epi64 (counter,0);
    sum +=  _mm_extract_epi64 (counter,1);

    sum += *randomsource;

    return (int) (sum>>32);
}


// Rabin-Karp Hashing
uint32_t hashRabinKarp(const uint64_t * useless, const uint32_t *  string, const size_t length) {
    int sum = 0;
    const int someprime = 31;
    for(size_t i = 0; i < length; ++i) {
        sum = someprime * sum +  string[i] ;
    }
    return sum;
}

// This is similar to Rabin-Karp, but we avoid the multiplication
//D. J. Bernstein, CDB–Constant Database, http://cr.yp.to/cdb.html
uint32_t hashBernstein(const uint64_t * useless,const uint32_t *  string, const size_t length) {
    int sum = 0;
    const int L = 3;
    for(size_t i = 0; i < length; ++i) {
        sum= (( sum<< L) + sum) ^ string[i] ;
    }
    return sum;
}

// Fowler-Noll-Vo hashing
// L. C. Noll, Fowler/Noll/Vo Hash, http://www.isthe.com/chongo/tech/comp/fnv/
uint32_t hashFNV1(const uint64_t * useless,const uint32_t *  string, const size_t length) {
    int sum = 0;
    const uint32_t someprime = 31;
    for(size_t i = 0; i < length; ++i) {
        sum = ( someprime * sum) ^ string[i] ;
    }
    return sum;
}

// Fowler-Noll-Vo hashing
// L. C. Noll, Fowler/Noll/Vo Hash, http://www.isthe.com/chongo/tech/comp/fnv/
uint32_t hashFNV1a(const uint64_t * useless ,const uint32_t *  string, const size_t length) {
    int sum = 0;
    const uint32_t someprime = 31;
    for(size_t i = 0; i < length; ++i) {
        sum = ( (string[i]) ^ sum) * someprime ;
    }
    return sum;
}

// M. V. Ramakrishna, J. Zobel, Performance in practice of string hashing functions,
//in: Proc. Int. Conf. on Database Systems for Advanced Applications, 1997.
uint32_t hashSAX(const uint64_t * useless,  const uint32_t *  string, const size_t length) {
    int sum = 0;
    const int L = 3;
    const int R = 5;
    for(size_t i = 0; i < length; ++i) {
        sum = sum ^((sum<<L)+(sum>>R) + string[i]) ;
    }
    return sum;
}


typedef unsigned long long ticks;

// Taken from stackoverflow (see http://stackoverflow.com/questions/3830883/cpu-cycle-count-based-profiling-in-c-c-linux-x86-64)
// Can give nonsensical results on multi-core AMD processors.
ticks rdtsc() {
    unsigned int lo, hi;
    asm volatile (
        "cpuid \n" /* serializing */
        "rdtsc"
        : "=a"(lo), "=d"(hi) /* outputs */
        : "a"(0)           /* inputs */
        : "%ebx", "%ecx");     /* clobbers*/
    return ((unsigned long long)lo) | (((unsigned long long)hi) << 32);
}
ticks startRDTSC (void) {
    return rdtsc();
}

ticks stopRDTSCP (void) {
    return rdtsc();
}
// start and stop are as recommended by 
// Gabriele Paoloni, How to Benchmark Code Execution Times on Intel® IA-32 and IA-64 Instruction Set Architectures
// September 2010
// http://edc.intel.com/Link.aspx?id=3954

static __inline__ ticks fancystartRDTSC (void) {
  unsigned cycles_low, cycles_high;
  asm volatile ("CPUID\n\t"
                "RDTSC\n\t"
                "mov %%edx, %0\n\t"
                "mov %%eax, %1\n\t": "=r" (cycles_high), "=r" (cycles_low)::
                "%rax", "%rbx", "%rcx", "%rdx");
  return ((ticks)cycles_high << 32) | cycles_low;
}

static __inline__ ticks fancystopRDTSCP (void) {
  unsigned cycles_low, cycles_high;
/// This should work fine on most machines, if the RDTSCP thing
/// fails for you, use the  rdtsc() call instead.
  asm volatile("RDTSCP\n\t"
               "mov %%edx, %0\n\t"
               "mov %%eax, %1\n\t"
               "CPUID\n\t": "=r" (cycles_high), "=r" (cycles_low):: "%rax",
               "%rbx", "%rcx", "%rdx");
  return ((ticks)cycles_high << 32) | cycles_low;
}


typedef uint32_t (*hashFunction)(const uint64_t *  ,const  uint32_t * , const size_t );

#define HowManyFunctions 6

hashFunction funcArr[HowManyFunctions] = {&hashMultilinear, &hashMultilinearSSE41, &hashRabinKarp, &hashFNV1, &hashFNV1a, &hashSAX};

 const char* functionnames[HowManyFunctions] = {"Multilinear   (strongly universal)",
 												"MultilinearSSE(strongly universal)", 
                                                 "RabinKarp                        ",
                                                 "FNV1                             ",
                                                 "FNV1a                            ",
                                                 "SAX                              "};

int main(int c, char ** arg) {
    const int N = 1024; // should be divisible by two!
    const int  SHORTTRIALS = 1000000;
    const int HowManyRepeats = 3;
    int i,k,j;
    int elapsed;
    hashFunction thisfunc;
    const char * functionname;
    ticks bef,aft;
    struct timeval start, finish;
    uint64_t randbuffer[N + 3];  
    uint32_t sumToFoolCompiler = 0;
    uint32_t intstring[N];// // could force 16-byte alignment with  __attribute__ ((aligned (16)));
    for (i = 0; i < N + 3; ++i) {
        randbuffer[i]=rand()| ((uint64_t)(rand())<<32);
    }
    for ( i = 0; i < N; ++i) {
        intstring[i] = rand();
    }
    printf("For documentation, see Strongly universal string hashing is fast at http://arxiv.org/abs/1202.4961");

    printf("Reporting the number of cycles per byte and the billions of bytes processed per second.\n");
    for(k = 0; k<HowManyRepeats; ++k) {
        printf("test #%d\n",k+1);
        for(i=0; i<HowManyFunctions; ++i) {
            thisfunc = funcArr[i];
            functionname = functionnames[i];
            gettimeofday( &start, 0);
            bef = startRDTSC();
            for(j=0; j < SHORTTRIALS; ++j)
                sumToFoolCompiler += thisfunc( &randbuffer[0],&intstring[0], N);
            aft = stopRDTSCP();
            gettimeofday( &finish, 0);  
            elapsed = ( 1000000*(finish.tv_sec - start.tv_sec) + (finish.tv_usec - start.tv_usec));
            printf("%s CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",functionname,
              (aft-bef)*1.0/(4.0*SHORTTRIALS*N),(4.0*SHORTTRIALS*N)/(1000.*elapsed));
        }
        printf("\n");
    }
    printf("# ignore this #%d\n",sumToFoolCompiler);

}


