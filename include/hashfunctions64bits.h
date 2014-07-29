

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
	/*
	 * We do not initialize ctx, and effectively ignore the random keys (rs), but we should
	 * really call:
	 *
	 * vmac_set_key(rs, &ctx)
	 *
	 * Doing so each time would almost surely handicap the speed of the function, so
	 * we ignore this.
	 */
	uint64_t tagl;// I think that this is useless but I am not sure, still needed due to API
	/*
	 * If the input is divisible by = VMAC_NHBYTES16 bytes, then we can call VHASH directly
	 * otherwise we have to do something messy due to alignment requirements in VHASH.
	 */
	size_t inputlengthinbytes = length * sizeof(uint64_t);
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
	}
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







/*******
 * Next macros are take directly from the original VMAC/VHASH implementation
 * by Ted Krovetz
 ******/

/* ----------------------------------------------------------------------- */
#if (__GNUC__ && (__x86_64__ || __amd64__))
/* ----------------------------------------------------------------------- */

#define ADD128(rh,rl,ih,il)                                               \
    asm ("addq %3, %1 \n\t"                                               \
         "adcq %2, %0"                                                    \
    : "+r"(rh),"+r"(rl)                                                   \
    : "r"(ih),"r"(il) : "cc");

#define MUL64(rh,rl,i1,i2)                                                \
    asm ("mulq %3" : "=a"(rl), "=d"(rh) : "a"(i1), "r"(i2) : "cc")

#define PMUL64 MUL64

#define GET_REVERSED_64(p)                                                \
    ({uint64_t x;                                                         \
     asm ("bswapq %0" : "=r" (x) : "0"(*(uint64_t *)(p))); x;})

/* ----------------------------------------------------------------------- */
#elif (__GNUC__ && __i386__)
/* ----------------------------------------------------------------------- */

#define GET_REVERSED_64(p)                                                \
    ({ uint64_t x;                                                        \
    uint32_t *tp = (uint32_t *)(p);                                       \
    asm  ("bswap %%edx\n\t"                                               \
          "bswap %%eax"                                                   \
    : "=A"(x)                                                             \
    : "a"(tp[1]), "d"(tp[0]));                                            \
    x; })

/* ----------------------------------------------------------------------- */
#elif (__GNUC__ && __ppc64__)
/* ----------------------------------------------------------------------- */

#define ADD128(rh,rl,ih,il)                                               \
    asm volatile (  "addc %1, %1, %3 \n\t"                                \
                    "adde %0, %0, %2"                                     \
    : "+r"(rh),"+r"(rl)                                                   \
    : "r"(ih),"r"(il));

#define MUL64(rh,rl,i1,i2)                                                \
{ uint64_t _i1 = (i1), _i2 = (i2);                                        \
    rl = _i1 * _i2;                                                       \
    asm volatile ("mulhdu %0, %1, %2" : "=r" (rh) : "r" (_i1), "r" (_i2));\
}

#define PMUL64 MUL64

#define GET_REVERSED_64(p)                                                \
    ({ uint32_t hi, lo, *_p = (uint32_t *)(p);                            \
       asm volatile ("lwbrx %0, %1, %2" : "=r"(lo) : "b%"(0), "r"(_p) );  \
       asm volatile ("lwbrx %0, %1, %2" : "=r"(hi) : "b%"(4), "r"(_p) );  \
       ((uint64_t)hi << 32) | (uint64_t)lo; } )

/* ----------------------------------------------------------------------- */
#elif (__GNUC__ && (__ppc__ || __PPC__))
/* ----------------------------------------------------------------------- */

#define GET_REVERSED_64(p)                                                \
    ({ uint32_t hi, lo, *_p = (uint32_t *)(p);                            \
       asm volatile ("lwbrx %0, %1, %2" : "=r"(lo) : "b%"(0), "r"(_p) );  \
       asm volatile ("lwbrx %0, %1, %2" : "=r"(hi) : "b%"(4), "r"(_p) );  \
       ((uint64_t)hi << 32) | (uint64_t)lo; } )

/* ----------------------------------------------------------------------- */
#elif (__GNUC__ && (__ARMEL__ || __ARM__))
/* ----------------------------------------------------------------------- */

#define bswap32(v)                                                        \
({ uint32_t tmp,out;                                                      \
    asm volatile(                                                         \
        "eor    %1, %2, %2, ror #16\n"                                    \
        "bic    %1, %1, #0x00ff0000\n"                                    \
        "mov    %0, %2, ror #8\n"                                         \
        "eor    %0, %0, %1, lsr #8"                                       \
    : "=r" (out), "=&r" (tmp)                                             \
    : "r" (v));                                                           \
    out;})

/* ----------------------------------------------------------------------- */
#elif _MSC_VER
/* ----------------------------------------------------------------------- */

#include <intrin.h>

#if (_M_IA64 || _M_X64) && \
    (!defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1000)
#define MUL64(rh,rl,i1,i2)   (rl) = _umul128(i1,i2,&(rh));
#pragma intrinsic(_umul128)
#define PMUL64 MUL64
#endif

/* MSVC uses add, adc in this version */
#define ADD128(rh,rl,ih,il)                                          \
    {   uint64_t _il = (il);                                         \
        (rl) += (_il);                                               \
        (rh) += (ih) + ((rl) < (_il));                               \
    }

#if _MSC_VER >= 1300
#define GET_REVERSED_64(p) _byteswap_uint64(*(uint64_t *)(p))
#pragma intrinsic(_byteswap_uint64)
#endif

#if _MSC_VER >= 1400 && \
    (!defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1000)
#define MUL32(i1,i2)    (__emulu((uint32_t)(i1),(uint32_t)(i2)))
#pragma intrinsic(__emulu)
#endif

/* ----------------------------------------------------------------------- */
#endif
/* ----------------------------------------------------------------------- */


/*******
 * End of bit taken from the VMAC/VHASH implementation by Ted Krovetz
 */


// like NH (from VHash) but without the pyramidal part nor any modulo reduction
uint64_t hashNH64(const void*  rs, const uint64_t *  string, const size_t length) {
	const uint64_t*  randomsource = (const uint64_t*) rs;
	uint64_t low = 0;
	uint64_t high = 0;
	size_t i = 0 ;
	for (; i + 7 < length; i+= 8) {
		uint64_t thigh, tlow;
	    MUL64(thigh,tlow,string[i]+randomsource[i],string[i+1]+randomsource[i+1]);
	    ADD128(high,low,thigh,tlow);
	    MUL64(thigh,tlow,string[i+2]+randomsource[i+2],string[i+3]+randomsource[i+3]);
	    ADD128(high,low,thigh,tlow);
	    MUL64(thigh,tlow,string[i+4]+randomsource[i+4],string[i+7]+randomsource[i+5]);
	    ADD128(high,low,thigh,tlow);
	    MUL64(thigh,tlow,string[i+6]+randomsource[i+6],string[i+7]+randomsource[i+7]);
	    ADD128(high,low,thigh,tlow);
	}
	for (; i + 1 < length; i+= 2) {
			uint64_t thigh, tlow;
		    MUL64(thigh,tlow,string[i]+randomsource[i],string[i+1]+randomsource[i+1]);
		    ADD128(high,low,thigh,tlow);
	}
	if(i < length) {
		uint64_t thigh, tlow;
	    MUL64(thigh,tlow,string[i]+randomsource[i],0+randomsource[i+1]);
	    ADD128(high,low,thigh,tlow);
	}
    // we use the prime 2^64 + 13
    return low + high * 13;// TODO: vhash actual computation is more complicated (and more expensive)
}

#endif /* HASHFUNCTIONS64BIT_H_ */
