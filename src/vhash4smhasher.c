#include <stdlib.h>
#include "vmac.h"
#include <string.h>

uint64_t vrandoms[140];
vmac_ctx_t ctx;

void init_vhash( uint32_t seed) { 
  srand(seed);
  for (int i=0; i < 140; ++i)
    vrandoms[i] = rand() | ( ((uint64_t) rand()) << 32);
  vmac_set_key((unsigned char *)vrandoms, &ctx);
}

uint64_t vhash4smhasher( const void *key, int lengthInBytes) {
  uint64_t tagl;// I think that this is useless but I am not sure, still needed due to API
  size_t inputlengthinbytes = lengthInBytes;
  return vhash((unsigned char *)key, inputlengthinbytes, &tagl, &ctx);
}



