
/**
 * Usage:
 *
 * $ make nhvsclnh.o
 *  $ /opt/intel/iaca/bin/iaca.sh -64 -mark 1 ./nhvsclnh.o
 * $ /opt/intel/iaca/bin/iaca.sh -64 -mark 2 ./nhvsclnh.o
 */


#include <stdint.h>
#include <stdlib.h>
#include <x86intrin.h>
#include </opt/intel/iaca/include/iacaMarks.h>

#define ADD128(rh,rl,ih,il)                                               \
    asm ("addq %3, %1 \n\t"                                               \
         "adcq %2, %0"                                                    \
    : "+r"(rh),"+r"(rl)                                                   \
    : "r"(ih),"r"(il) : "cc");

#define MUL64(rh,rl,i1,i2)                                                \
    asm ("mulq %3" : "=a"(rl), "=d"(rh) : "a"(i1), "r"(i2) : "cc")

void testNH(uint64_t * rh, uint64_t * rl, uint64_t * m, uint64_t * k ) {

	IACA_START;
	uint64_t u = m[0]+k[0];
	uint64_t v = m[1]+k[1];
	__asm__ ("mulq %[v]\n"
             "addq %%rax,  %[rl]\n"
              "adcq %%rdx,  %[rh]\n"
             :  [rh] "+m" (*rh), [rl] "+m" (*rl)  : [u] "a" (u), [v] "r" (v)  :"rdx","cc" );
		IACA_END;
}

void testCLNH(__m128i * acc, uint64_t * m, uint64_t * k ) {
	IACA_START;
	const __m128i temp1 = _mm_load_si128((__m128i *) k);
	const __m128i temp2 = _mm_lddqu_si128((__m128i *) m);
	const __m128i add1 = _mm_xor_si128(temp1, temp2);
	const __m128i clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
	*acc = _mm_xor_si128(clprod1, *acc);
	IACA_END;
}
