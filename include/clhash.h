
#ifndef CLHASHCONVENIENCE_H_
#define CLHASHCONVENIENCE_H_



/////////
// what follows are convenience functions
// call init_clhash once with a 32-bit key
// then call clhash to hash strings.
/////////
#include "mersenne.h"
#include "clmulhierarchical64bits.h"

static uint64_t randomkey[RANDOM_64BITWORDS_NEEDED_FOR_CLHASH];

void init_clhash( uint32_t seed) {
    ZRandom zr;
    initZRandom(&zr,seed);
    for (int i=0; i < RANDOM_64BITWORDS_NEEDED_FOR_CLHASH; ++i)
        randomkey[i] = getValue(&zr) | ( ((uint64_t) getValue(&zr)) << 32);
}

uint64_t clhash( const void *key, int len) {
    return CLHASHbyte(randomkey,(const char *)key,len);
}



#endif /* CLHASHCONVENIENCE_H_ */
