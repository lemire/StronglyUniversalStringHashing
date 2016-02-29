

#ifndef HASHFUNCTIONS64BIT_H_
#define HASHFUNCTIONS64BIT_H_

// first pointer is a random source,
// next pointer is the data.
// outputs a hash value
typedef uint64_t (*hashFunction64)(const void *  ,const  uint64_t * , const size_t );


#include "PMP/PMP_C_wrapper.h"

// PMP64 hash function (ignores seed, but good to benchmark running speed)
uint64_t hashPMP64(const void*  rs, const uint64_t *  string, const size_t length) {
    (void) rs;
    return pmp64_hash((const unsigned char *) string, length * sizeof(uint64_t));
}


#include "City/City.h"

// Google hash function
uint64_t hashCity(const void*  rs, const uint64_t *  string, const size_t length) {
    return CityHash64WithSeed((const char *) string, length * sizeof(uint64_t),*(uint64_t*)rs);
}

#include "SipHash/siphash24.h"
// SipHash
uint64_t hashSipHash(const void*  rs, const uint64_t *  string, const size_t length) {
    uint64_t answer;
    siphash((uint8_t * ) &answer, (const uint8_t * )string, length *sizeof(uint64_t), (const uint8_t *)rs );
    return answer;
}



#include "VHASH/vmac.h"
#include <string.h>
// to simulate the speed of VHASH. Not thread safe.
// This is not correct if the input is not divisible by 16 bytes.
uint64_t hashVHASH64(const void*  rs, const uint64_t *  string, const size_t length) {
    static vmac_ctx_t ctx; // we need a ctx struct
    /*
     * At least once, we need the following call:
     *
     * vmac_set_key(rs, &ctx)
     *
     * Doing so each time would almost surely handicap the speed of the function, so
     * we only do it if our heuristic suggests that it is necessary.
     *
     * This heuristic can be defeated: (1) we could end up calling vmac_set_key
     * unnecessarily (2) we could fail to call it even though it is necessary.
     * Neither of these problems will occur in our benchmark.
     *
     *
     */
    static const void *initialized_rs = 0;
    static uint64 first_value = 0;

    //  initialize, hoping a change of randoms will change the first value or the buffer location
    if (rs != initialized_rs || (first_value != *((uint64*) rs))) {
        vmac_set_key((unsigned char *) rs, &ctx);
        initialized_rs = rs;
        first_value = *((uint64*) rs);
    }

    uint64_t tagl;// I think that this is useless but I am not sure, still needed due to API
    /*
     * If the input is divisible by = VMAC_NHBYTES16 bytes, then we can call VHASH directly
     * otherwise we have to do something messy due to alignment requirements in VHASH.
     */
    size_t inputlengthinbytes = length * sizeof(uint64_t);
    return vhash((unsigned char *)string, inputlengthinbytes, &tagl, &ctx);
}

// to simulate the speed of VHASH. Not thread safe.
// This is  correct even the input is not divisible by 16 bytes, but slower than hashVHASH64.
uint64_t saferhashVHASH64(const void*  rs, const uint64_t *  string, const size_t length) {
    static vmac_ctx_t ctx; // we need a ctx struct
    /*
     * At least once, we need the following call:
     *
     * vmac_set_key(rs, &ctx)
     *
     * Doing so each time would almost surely handicap the speed of the function, so
     * we only do it if our heuristic suggests that it is necessary.
     *
     * This heuristic can be defeated: (1) we could end up calling vmac_set_key
     * unnecessarily (2) we could fail to call it even though it is necessary.
     * Neither of these problems will occur in our benchmark.
     *
     *
     */
    static const void *initialized_rs = 0;
    static uint64 first_value = 0;

    //  initialize, hoping a change of randoms will change the first value or the buffer location
    if (rs != initialized_rs || (first_value != *((uint64*) rs))) {
        vmac_set_key((unsigned char *) rs, &ctx);
        initialized_rs = rs;
        first_value = *((uint64*) rs);
    }
    uint64_t tagl;// I think that this is useless but I am not sure, still needed due to API
    /*
     * If the input is divisible by = VMAC_NHBYTES16 bytes, then we can call VHASH directly
     * otherwise we have to do something messy due to alignment requirements in VHASH.
     */
    size_t inputlengthinbytes = length * sizeof(uint64_t);
    return vhash((unsigned char *)string, inputlengthinbytes, &tagl, &ctx);
    if((inputlengthinbytes % VMAC_NHBYTES) == 0) {
        return vhash((unsigned char *)string, inputlengthinbytes, &tagl, &ctx);
    } else {
        size_t roundedlength = inputlengthinbytes/VMAC_NHBYTES*VMAC_NHBYTES;
        vhash_update((unsigned char*)string, roundedlength, &ctx);
        if(roundedlength > 0) vhash_update((unsigned char*)string, roundedlength, &ctx);
        static unsigned char lastBlock[VMAC_NHBYTES + 16];
        unsigned char *alignedptr = (unsigned char*)(((uintptr_t)lastBlock+15) & ~ (uintptr_t)0x0F);
        size_t remaining = inputlengthinbytes - roundedlength;
        memcpy(alignedptr, (unsigned char*)string + roundedlength, remaining);
        memset(alignedptr + remaining, 0, VMAC_NHBYTES - remaining );
        return vhash(alignedptr, remaining, &tagl, &ctx);
    }
}




#endif /* HASHFUNCTIONS64BIT_H_ */
