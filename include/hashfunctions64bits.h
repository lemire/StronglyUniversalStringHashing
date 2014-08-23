

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



#include "vmac.h"
#include <string.h>
// to simulate the speed of VHASH. Not thread safe.
uint64_t hashVHASH64(const void*  rs, const uint64_t *  string, const size_t length) {
	static vmac_ctx_t ctx; // we need a ctx struct
        static void *initialized_rs = 0; // temp temp
        static uint64 first_value = 0; // temp temp
	/*
	 * We do not initialize ctx, and effectively ignore the random keys (rs), but we should
	 * really call:
	 *
	 * vmac_set_key(rs, &ctx)
	 *
	 * Doing so each time would almost surely handicap the speed of the function, so
	 * we ignore this.
	 */
      
        // temp temp, initialize, hoping a change of randoms will change the first value or the buffer location 
        if (rs != initialized_rs || (first_value != *( (uint64*)rs) )) {
          vmac_set_key((unsigned char *)rs, &ctx);
          initialized_rs = rs;
          first_value = *( (uint64*) rs);
        }
          
	uint64_t tagl;// I think that this is useless but I am not sure, still needed due to API
	/*
	 * If the input is divisible by = VMAC_NHBYTES16 bytes, then we can call VHASH directly
	 * otherwise we have to do something messy due to alignment requirements in VHASH.
	 */
	size_t inputlengthinbytes = length * sizeof(uint64_t);
	return vhash((unsigned char *)string, inputlengthinbytes, &tagl, &ctx);
	// code below might not be necessary?
/*
	if((inputlengthinbytes % VMAC_NHBYTES) == 0) {
		return vhash((unsigned char *)string, inputlengthinbytes, &tagl, &ctx);
	} else {
		size_t roundedlength = inputlengthinbytes/VMAC_NHBYTES*VMAC_NHBYTES;
		if(roundedlength > 0) vhash_update((unsigned char*)string, roundedlength, &ctx);
		unsigned char lastBlock[VMAC_NHBYTES + 16];
		unsigned char *alignedptr = (unsigned char*)(((uintptr_t)lastBlock+15) & ~ (uintptr_t)0x0F);
		size_t remaining = inputlengthinbytes - roundedlength;
		assert(remaining + roundedlength == inputlengthinbytes);
		memcpy(alignedptr, (unsigned char*)string + roundedlength, remaining);
		memset(alignedptr + remaining, 0, VMAC_NHBYTES - remaining );
		return vhash(alignedptr, remaining, &tagl, &ctx);
	}*/
}





#endif /* HASHFUNCTIONS64BIT_H_ */
