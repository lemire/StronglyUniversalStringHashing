#include <string.h>

#include "clmulhierarchical64bits.h"
#include "mersenne.h"

// input is arbitrary length (may have to pad)


static uint64_t randoms[RANDOM_64BITWORDS_NEEDED_FOR_CLHASH];

void init_cl3264( uint32_t seed) {
    ZRandom zr;
    initZRandom(&zr,seed);
    for (int i=0; i < RANDOM_64BITWORDS_NEEDED_FOR_CLHASH; ++i)
        randoms[i] = getValue(&zr) | ( ((uint64_t) getValue(&zr)) << 32);
}

// 64 bit version, for comparison vs cityhash with smhasher
uint64_t cl64( const void *key, int len) {
    return CLHASHbyte(randoms,key,len);
}


// define a 32-bit version of CLMUL for quality tests.  Just use low-order bits.
uint32_t cl32( const void *key, int len) {
    return (uint32_t) cl64(key,len); // low bits
}
