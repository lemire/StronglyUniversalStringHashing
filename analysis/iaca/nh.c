/**
 * Usage:
 *
 * $ make nhvsclnh.o
 * $ /opt/intel/iaca/bin/iaca.sh -64 -mark 1 ./nhvsclnh.o
 * $ /opt/intel/iaca/bin/iaca.sh -64 -mark 2 ./nhvsclnh.o
 */

#include <stdint.h>
#include <stdlib.h>
#include <x86intrin.h>
#ifdef IACA
#include </opt/intel/iaca/include/iacaMarks.h>
#else
#define IACA_START
#define IACA_END
#endif


void testNH(uint64_t * rh, uint64_t * rl, uint64_t * m, uint64_t * k ) {
    IACA_START;
    uint64_t u = m[0]+k[0];
    uint64_t v = m[1]+k[1];
    __asm__ ("mulq %[v]\n"
             "addq %%rax,  %[rl]\n"
             "adcq %%rdx,  %[rh]\n"
             :  [rh] "+m" (*rh), [rl] "+m" (*rl)  : [u] "a" (u), [v] "r" (v)  :"rdx","cc" );
    u = m[2]+k[2];
    v = m[3]+k[3];
    __asm__ ("mulq %[v]\n"
             "addq %%rax,  %[rl]\n"
             "adcq %%rdx,  %[rh]\n"
             :  [rh] "+m" (*rh), [rl] "+m" (*rl)  : [u] "a" (u), [v] "r" (v)  :"rdx","cc" );
    u = m[4]+k[4];
    v = m[5]+k[5];
    __asm__ ("mulq %[v]\n"
             "addq %%rax,  %[rl]\n"
             "adcq %%rdx,  %[rh]\n"
             :  [rh] "+m" (*rh), [rl] "+m" (*rl)  : [u] "a" (u), [v] "r" (v)  :"rdx","cc" );
    u = m[6]+k[6];
    v = m[7]+k[7];
    __asm__ ("mulq %[v]\n"
             "addq %%rax,  %[rl]\n"
             "adcq %%rdx,  %[rh]\n"
             :  [rh] "+m" (*rh), [rl] "+m" (*rl)  : [u] "a" (u), [v] "r" (v)  :"rdx","cc" );

    IACA_END;
}

void testCLNH(__m128i * acc, uint64_t * m, uint64_t * k ) {
    IACA_START;
    __m128i temp1, temp2, add1, clprod1;
    temp1 = _mm_load_si128((__m128i *) k);
    temp2 = _mm_lddqu_si128((__m128i *) m);
    add1 = _mm_xor_si128(temp1, temp2);
    clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
    *acc = _mm_xor_si128(clprod1, *acc);
    temp1 = _mm_load_si128((__m128i *) k + 1);
    temp2 = _mm_lddqu_si128((__m128i *) m + 1);
    add1 = _mm_xor_si128(temp1, temp2);
    clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
    *acc = _mm_xor_si128(clprod1, *acc);
    temp1 = _mm_load_si128((__m128i *) k + 2);
    temp2 = _mm_lddqu_si128((__m128i *) m + 2);
    add1 = _mm_xor_si128(temp1, temp2);
    clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
    *acc = _mm_xor_si128(clprod1, *acc);
    temp1 = _mm_load_si128((__m128i *) k + 3);
    temp2 = _mm_lddqu_si128((__m128i *) m + 3);
    add1 = _mm_xor_si128(temp1, temp2);
    clprod1 = _mm_clmulepi64_si128(add1, add1, 0x10);
    *acc = _mm_xor_si128(clprod1, *acc);
    IACA_END;
}
