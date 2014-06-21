

#ifndef HASHFUNCTIONS64BIT_H_
#define HASHFUNCTIONS64BIT_H_

#include "City.h"
// first pointer is a random source,
// next pointer is the data.
// outputs a hash value
typedef uint64_t (*hashFunction64)(const void *  ,const  uint64_t * , const size_t );


// Google hash function
uint64_t hashCity(const void*  rs, const uint64_t *  string, const size_t length) {
	return CityHash64WithSeed((const char *) string, length * sizeof(uint64_t),*(uint64_t*)rs);
}

// like MHH, this is essentially multilinear with 64bit multiplication
// summed up over a 128-bit counter
uint64_t hashMMH_NonPyramidal(const void*  rs, const uint64_t *  string, const size_t length) {
	const uint64_t*  randomsource = (const uint64_t*) rs;

    uint64_t low = 0;
    uint64_t high = 0;
    size_t i = 0;
    for(; i<length/8*8; i+=8) {
    __asm__ (
    "movq (%[u]),%%rax\n"
    "mulq (%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 8(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 8(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 16(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 16(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 24(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 24(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 32(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 32(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 40(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 40(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 48(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 48(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "movq 56(%[u]),%%rax\n"
    "adcq %%rdx,  %[rh]\n"
    "mulq 56(%[v])\n"
    "addq %%rax,  %[rl]\n"
    "adcq %%rdx,  %[rh]\n"
                 :  [rh] "+r" (high) , [rl] "+r" (low)  : [u] "r" (randomsource+i), [v] "r" (string+i)  :"rdx","rax","memory","cc");
    }

    for(; i<length; ++i) {
        __asm__ ("mulq %[v]\n"
                 "addq %%rax,  %[rl]\n"
                 "adcq %%rdx,  %[rh]\n"
                 :  [rh] "+r" (high), [rl] "+r" (low)  : [u] "a" (randomsource[i]), [v] "r" (string[i])  :"rdx","cc");
    }
    // we use the prime 2^64 + 13
    return low + high * 13;// TODO: actual computation is more complicated (and more expensive) so we slightly help MMH here
}


#endif /* HASHFUNCTIONS64BIT_H_ */
