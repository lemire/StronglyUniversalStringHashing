#include <string.h>

#include "clmulhierarchical64bits.h"
#include "mersenne.h"

// input is arbitrary length (may have to pad)

uint64_t randoms[150];

void init_cl3264( uint32_t seed) { 
  ZRandom zr;
  initZRandom(&zr,seed);
  for (int i=0; i < 150; ++i)
    randoms[i] = getValue(&zr) | ( ((uint64_t) getValue(&zr)) << 32);
}

// 64 bit version, for comparison vs cityhash with smhasher
uint64_t cl64( const void *key, int len) {
	return CLHASHbyte(randoms,key,len);
/*  void *keypadded;
  int newlen = len;
  if ( (len & 15) == 0)
    keypadded = (void *)key;
  else {  // this is going to hurt, we could do better
    newlen = (len+15) & ~15;
    keypadded = malloc( newlen);
    memcpy(keypadded, key, len);
    for (int l = len; l < newlen; ++l)
      // pad with ones; one of the smhasher sanity tests was to verify that
      // user padding with zeros always changes the hash value, and
      // we also pad with zeros, the sanity check fails.
      * ( (char *) keypadded+l) = 1;
  }
  uint64_t ans = CLHASH( randoms, keypadded,newlen/8);  // newlen in 64-bits
  if (keypadded != key)  // we had malloc'ed
    free(keypadded);
  return ans;*/
}


// define a 32-bit version of CLMUL for quality tests.  Just use low-order bits.
uint32_t cl32( const void *key, int len) {
  return (uint32_t) cl64(key,len); // low bits
}
