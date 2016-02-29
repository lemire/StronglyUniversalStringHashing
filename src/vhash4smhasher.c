#include <stdlib.h>
#include <string.h>

#include "mersenne.h"
#include "VHASH/vmac.h"


uint64_t vrandoms[140];
vmac_ctx_t ctx;

void init_vhash( uint32_t seed) {
    ZRandom zr;
    initZRandom(&zr,seed);
    for (int i=0; i < 140; ++i)
        vrandoms[i] = getValue(&zr) | ( ((uint64_t) getValue(&zr)) << 32);
    vmac_set_key((unsigned char *)vrandoms, &ctx);
}

uint64_t vhash4smhasher( const void *key, int lengthInBytes) {
    uint64_t tagl;// I think that this is useless but I am not sure, still needed due to API
    size_t inputlengthinbytes = lengthInBytes;
    if((inputlengthinbytes % VMAC_NHBYTES) == 0) {
        return vhash((unsigned char *)key, inputlengthinbytes, &tagl, &ctx);
    }
    size_t roundedlength = inputlengthinbytes/VMAC_NHBYTES*VMAC_NHBYTES;
    if(roundedlength > 0) vhash_update((unsigned char*)key, roundedlength, &ctx);
    unsigned char lastBlock[VMAC_NHBYTES + 16];
    unsigned char *alignedptr = (unsigned char*)(((uintptr_t)lastBlock+15) & ~ (uintptr_t)0x0F);
    size_t remaining = inputlengthinbytes - roundedlength;
    memcpy(alignedptr, (unsigned char*)key + roundedlength, remaining);
    memset(alignedptr + remaining, 0, VMAC_NHBYTES - remaining );
    return vhash(alignedptr, remaining, &tagl, &ctx);
}

