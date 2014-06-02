
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
 
//
// this is strongly universal. Condition: randomsource must be at least as long as 
// the string length  + 3.
//  
// Reference: Owen Kaser and Daniel Lemire, Strongly universal string hashing is fast, Computer Journal 
// http://arxiv.org/abs/1202.4961
#if defined(__GNUC__) && !( defined(__clang__) || defined(__INTEL_COMPILER  ))
__attribute__((optimize("no-tree-vectorize")))
// GCC has buggy SSE2 code generation in some cases
// Thanks to Nathan Kurz for noticing that GCC 4.7 requires no-tree-vectorize to produce correct results.
#endif
uint32_t hashMultilinear(const uint64_t *  randomsource, const uint32_t *  string, const size_t length) {
    const uint32_t * const endstring = string + length;
    uint64_t sum = *(randomsource++);
    for(; string!= endstring; ++randomsource,++string ) {
        sum+= (*randomsource *  (uint64_t)(*string)) ;
    }
    sum += *randomsource;
    return (int) (sum>>32);
}

#if defined(__GNUC__) && !( defined(__clang__) || defined(__INTEL_COMPILER  ))
__attribute__((optimize("no-tree-vectorize")))
// GCC has buggy SSE2 code generation in some cases
#endif
uint32_t hashMultilinear2by2(const uint64_t *  randomsource, const uint32_t *  string, const size_t length) {
    assert ( length / 2 * 2 == length );// length is even
    const uint32_t * const endstring = string + length;
    uint64_t sum = *(randomsource++);
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        sum+= (*randomsource *  (uint64_t)(*string)) + (*(randomsource+1) *  (uint64_t)(*(string+1)))  ;
    }
    sum += *randomsource;
    return (int) (sum>>32);
}

#if defined(__GNUC__) && !( defined(__clang__) || defined(__INTEL_COMPILER  ))
__attribute__((optimize("no-tree-vectorize")))
// GCC has buggy SSE2 code generation in some cases
#endif
uint32_t hashMultilinearhalf(const uint64_t *  randomsource, const uint32_t *  string, const size_t length) {
    assert ( length / 2 * 2 == length );// length is even
    const uint32_t * const endstring = string + length;
    uint64_t sum = *(randomsource++);
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        sum+= (*randomsource +  (uint64_t)(*string)) * (*(randomsource+1) +  (uint64_t)(*(string+1)))  ;
    }
    sum += *randomsource;
    return (int) (sum>>32);
}

#if defined(__GNUC__) && !( defined(__clang__) || defined(__INTEL_COMPILER  ))
__attribute__((optimize("no-tree-vectorize")))
// GCC has buggy SSE2 code generation in some cases
#endif
uint32_t hashMultilineardouble(const uint64_t *  randomsource, const uint32_t *  string, const size_t length) {
    assert ( length / 2 * 2 == length );// length is pair
    const uint32_t * const endstring = string + length;
    uint64_t sum = *(randomsource++);
    uint64_t s2 = 0;
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        sum += *randomsource *  (uint64_t)(*string);
        s2 +=  *(randomsource+1) *  (uint64_t)(*(string+1))  ;
    }
    sum += *randomsource + s2;
    return (int) (sum>>32);
}



// used by pyramidal_Multilinear below
// code not thoroughly checked
void __hashMulti(const uint64_t *  randomsource, const uint32_t *  string, uint32_t *  output, int length, int blocksize) {
    int bufferpos = 0;
    int offset = length - length/blocksize*blocksize;
    for(; bufferpos<length/blocksize;++bufferpos) {
          output[bufferpos] = hashMultilinear(randomsource, string+bufferpos*blocksize, blocksize);
    } 
    if(offset>0) {
        output[length/blocksize] = hashMultilinear(randomsource, string+length/blocksize*blocksize, 
        offset);
    }
} 


// this function is 4/2**32 almost universal on 32 bits
// it uses at most 8KB of random bits (for strings having 32-bit lengths)
// code not thoroughly checked
uint32_t pyramidal_Multilinear(const uint64_t *  randomsource, const uint32_t *  string, const size_t len) {
	size_t length = len;
    int blocksize = 256;  
    int newlength = ((length+blocksize-1)/blocksize);
    uint32_t * array = malloc(newlength * sizeof(uint32_t));
    __hashMulti(randomsource, string, array,  length, blocksize);
    randomsource+=blocksize+1;
    length = newlength;
    while(length>1) {
        newlength = ((length+blocksize-1)/blocksize);
        uint32_t * array2 = malloc(newlength * sizeof(uint32_t));
        __hashMulti(randomsource, array, array2,  length, blocksize);
        randomsource+=blocksize+1;
        free(array);
        array = array2;
        length = newlength;
        if(length == 1) {
        }
    }
    uint32_t answer = array[0];
    free(array);
    return answer;
}


// Rabin-Karp Hashing
uint32_t hashRabinKarp(const uint64_t * useless, const uint32_t *  string, const size_t length) {
    (void) (useless);
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
    (void) (useless);
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
    (void) (useless);
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
    (void) (useless);
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
    (void) (useless);
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
typedef uint64_t (*hashFunction64)(const uint64_t *  ,const  uint64_t * , const size_t );


#ifdef __PCLMUL__

#include "clmul.h"

#define HowManyFunctions 11
#define HowManyFunctions64 4

hashFunction64 funcArr64[HowManyFunctions64] = {&hashGaloisFieldfast64,&hashGaloisFieldfast64half,&hashGaloisFieldfast64halfunrolled,&referenceproduct};

hashFunction funcArr[HowManyFunctions] = {&hashGaloisFieldMultilinear,
 &hashGaloisFieldMultilinearHalfMultiplications, &hashMultilinear,&hashMultilinear2by2 ,
 &hashMultilinearhalf, &hashMultilineardouble,
 &hashRabinKarp, &hashFNV1, &hashFNV1a, &hashSAX,&pyramidal_Multilinear};

const char* functionnames64[HowManyFunctions64] = { 
                                                "GFMultilinear (64-bit proto)        ",
                                                "GFMultilinear (64-bit half)         ",
                                                "GFMultilinear (64-bit half, unrol)  ",
                                                "Reference (Like MHH)                "};
                                                
                                                const char* functionnames[HowManyFunctions] = {
                                                "GFMultilinear   (strongly universal)",
                                                "GFMultilinearhalf   (str. universal)",
                                                "Multilinear     (strongly universal)",
                                                "Multilinear2x2  (strongly universal)",
                                                "Multilinearhalf (strongly universal)",
                                                "Multilineardouble (strongly u.)     ",
                                                "RabinKarp                           ",
                                                "FNV1                                ",
                                                "FNV1a                               ",
                                                "SAX                                 ",
                                                "Pyramidal multilinear (a. univ.)    "};
#else



#define HowManyFunctions 9


hashFunction funcArr[HowManyFunctions] = {&hashMultilinear,&hashMultilinear2by2 ,
&hashMultilinearhalf, &hashMultilineardouble,
 &hashRabinKarp, &hashFNV1, &hashFNV1a, &hashSAX,&pyramidal_Multilinear};

const char* functionnames[HowManyFunctions] = {"Multilinear  (strongly universal)",
                                                "Multilinear2x2  (strongly universal)",
                                                "Multilinearhalf (strongly universal)",
                                                "Multilineardouble (strongly u.)     ",
                                                "RabinKarp                           ",
                                                "FNV1                                ",
                                                "FNV1a                               ",
                                                "SAX                                 ",
                                                "Pyramidal multilinear (a. univ.)    "};


#endif

int main(int c, char ** arg) {
    (void) (c);
    (void) (arg);
    const int N = 1024; // should be divisible by two!
    const int  SHORTTRIALS = 1000000;
    const int HowManyRepeats = 3;
    int i,k,j;
    int elapsed;
    hashFunction thisfunc;
    
    const char * functionname;
    ticks bef,aft;
    struct timeval start, finish;
    uint64_t randbuffer[N + 3]  __attribute__ ((aligned (16)));  
    uint32_t sumToFoolCompiler = 0;
    uint32_t intstring[N]  __attribute__ ((aligned (16)));// // could force 16-byte alignment with  __attribute__ ((aligned (16)));
    for (i = 0; i < N + 3; ++i) {
        randbuffer[i]=rand()| ((uint64_t)(rand())<<32);
    }
    for ( i = 0; i < N; ++i) {
        intstring[i] = rand();
    }
    printf("For documentation, see Strongly universal string hashing is fast at http://arxiv.org/abs/1202.4961 \n");

    printf("Reporting the number of cycles per byte and the billions of bytes processed per second.\n");
    for(k = 0; k<HowManyRepeats; ++k) {
        printf("test #%d (64-bit hash values)\n",k+1);
#ifdef __PCLMUL__
        hashFunction64 thisfunc64;
        for(i=0; i<HowManyFunctions64; ++i) {
            sumToFoolCompiler = 0;
            thisfunc64 = funcArr64[i];
            functionname = functionnames64[i];
            gettimeofday( &start, 0);
            bef = startRDTSC();
            assert(N/2*2==N);
            for(j=0; j < SHORTTRIALS; ++j)
                sumToFoolCompiler += thisfunc64( &randbuffer[0],(uint64_t *)&intstring[0], N/2);
            aft = stopRDTSCP();
            gettimeofday( &finish, 0);  
            elapsed = ( 1000000*(finish.tv_sec - start.tv_sec) + (finish.tv_usec - start.tv_usec));
            printf("%s CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",functionname,
              (aft-bef)*1.0/(4.0*SHORTTRIALS*N),(4.0*SHORTTRIALS*N)/(1000.*elapsed));
                  printf("# ignore this #%d\n",sumToFoolCompiler);

        }
        printf("\n");
#endif 
        printf("test #%d (32-bit hash values)\n",k+1);
        for(i=0; i<HowManyFunctions; ++i) {
            sumToFoolCompiler = 0;
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
                  printf("# ignore this #%d\n",sumToFoolCompiler);

        }
        printf("\n");
    }
    printf("# ignore this #%d\n",sumToFoolCompiler);

}


