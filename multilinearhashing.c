
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



//
// this is strongly universal. Condition: randomsource must be at least as long as 
// the string length  + 3.
//  
// Reference: Owen Kaser and Daniel Lemire, Strongly universal string hashing is fast, Computer Journal 
// http://arxiv.org/abs/1202.4961
__attribute__ ((__target__ ("no-sse2"))) // GCC has buggy SSE2 code generation
uint32_t hashMultilinear(const uint64_t *  randomsource, const uint32_t *  string, const size_t length) {
    uint64_t sum = randomsource[0];
    for(size_t i = 0; i < length; ++i) {
        sum += (randomsource[i+2] *  (uint64_t)(string[i])) ;
    }
    sum += randomsource[length+1]; // always assume last "character" is a one
    return (uint32_t) (sum>>32);
}

// Rabin-Karp Hashing
uint32_t hashRabinKarp(const uint64_t * , const uint32_t *  string, const size_t length) {
    int sum = 0;
    const int someprime = 31;
    for(size_t i = 0; i < length; ++i) {
        sum = someprime * sum +  string[i] ;
    }
    return sum;
}

// This is similar to Rabin-Karp, but we avoid the multiplication
//D. J. Bernstein, CDB–Constant Database, http://cr.yp.to/cdb.html
uint32_t hashBernstein(const uint64_t * ,const uint32_t *  string, const size_t length) {
    int sum = 0;
    const int L = 3;
    for(size_t i = 0; i < length; ++i) {
        sum= (( sum<< L) + sum) ^ string[i] ;
    }
    return sum;
}

// Fowler-Noll-Vo hashing
// L. C. Noll, Fowler/Noll/Vo Hash, http://www.isthe.com/chongo/tech/comp/fnv/
uint32_t hashFNV1(const uint64_t * ,const uint32_t *  string, const size_t length) {
    int sum = 0;
    const uint32_t someprime = 31;
    for(size_t i = 0; i < length; ++i) {
        sum = ( someprime * sum) ^ string[i] ;
    }
    return sum;
}

// Fowler-Noll-Vo hashing
// L. C. Noll, Fowler/Noll/Vo Hash, http://www.isthe.com/chongo/tech/comp/fnv/
uint32_t hashFNV1a(const uint64_t *  ,const uint32_t *  string, const size_t length) {
    int sum = 0;
    const uint32_t someprime = 31;
    for(size_t i = 0; i < length; ++i) {
        sum = ( (string[i]) ^ sum) * someprime ;
    }
    return sum;
}

// M. V. Ramakrishna, J. Zobel, Performance in practice of string hashing functions,
//in: Proc. Int. Conf. on Database Systems for Advanced Applications, 1997.
uint32_t hashSAX(const uint64_t * ,  const uint32_t *  string, const size_t length) {
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


typedef uint32_t (*hashFunction)(const uint64_t *  ,const  uint32_t * , const size_t );

#define HowManyFunctions 5

hashFunction funcArr[HowManyFunctions] = {&hashMultilinear, &hashRabinKarp, &hashFNV1, &hashFNV1a, &hashSAX};

 const char* functionnames[HowManyFunctions] = {"Multilinear  (strongly universal)", 
 												"RabinKarp                        ",
 												"FNV1                             ",
 												"FNV1a                            ",
 												"SAX                              "};

int main(int , char **) {
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
    for (i=0; i < N + 2; ++i) {
        randbuffer[i]=rand()| ((uint64_t)(rand())<<32);
    }
    for (i=0; i < N + 1; ++i) {
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
            printf("%s CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",functionname,(aft-bef)*1.0/(4.0*SHORTTRIALS*N),(4.0*SHORTTRIALS*N)/(1000.*elapsed));
        }
        printf("\n");
    }
    printf("# ignore this #%d\n",sumToFoolCompiler);

}


