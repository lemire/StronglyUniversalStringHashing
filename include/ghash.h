/*
 * GHASH implementation based "Intel Carry-Less Multiplication Instruction and
 * its Usage for Computing the GCM Mode" by Shay Gueron, Michael E. Kounavis.
 */

#ifndef GHASH_H_
#define GHASH_H_

#include "clmul.h"

// Gueron and Kounavis fig. 5
__m128i gfmul_fig5(__m128i a, __m128i b) {
    __m128i tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
    tmp3 = _mm_clmulepi64_si128(a, b, 0x00);
    tmp4 = _mm_clmulepi64_si128(a, b, 0x10);
    tmp5 = _mm_clmulepi64_si128(a, b, 0x01);
    tmp6 = _mm_clmulepi64_si128(a, b, 0x11);

    tmp4 = _mm_xor_si128(tmp4, tmp5);
    tmp5 = _mm_slli_si128(tmp4, 8);
    tmp4 = _mm_srli_si128(tmp4, 8);
    tmp3 = _mm_xor_si128(tmp3, tmp5);
    tmp6 = _mm_xor_si128(tmp6, tmp4);

    tmp7 = _mm_srli_epi32(tmp3, 31);
    tmp8 = _mm_srli_epi32(tmp6, 31);
    tmp3 = _mm_slli_epi32(tmp3, 1);
    tmp6 = _mm_slli_epi32(tmp6, 1);

    tmp9 = _mm_srli_si128(tmp7, 12);
    tmp8 = _mm_slli_si128(tmp8, 4);
    tmp7 = _mm_slli_si128(tmp7, 4);
    tmp3 = _mm_or_si128(tmp3, tmp7);
    tmp6 = _mm_or_si128(tmp6, tmp8);
    tmp6 = _mm_or_si128(tmp6, tmp9);

    tmp7 = _mm_slli_epi32(tmp3, 31);
    tmp8 = _mm_slli_epi32(tmp3, 30);
    tmp9 = _mm_slli_epi32(tmp3, 25);

    tmp7 = _mm_xor_si128(tmp7, tmp8);
    tmp7 = _mm_xor_si128(tmp7, tmp9);
    tmp8 = _mm_srli_si128(tmp7, 4);
    tmp7 = _mm_slli_si128(tmp7, 12);
    tmp3 = _mm_xor_si128(tmp3, tmp7);

    tmp2 = _mm_srli_epi32(tmp3, 1);
    tmp4 = _mm_srli_epi32(tmp3, 2);
    tmp5 = _mm_srli_epi32(tmp3, 7);
    tmp2 = _mm_xor_si128(tmp2, tmp4);
    tmp2 = _mm_xor_si128(tmp2, tmp5);
    tmp2 = _mm_xor_si128(tmp2, tmp8);
    tmp3 = _mm_xor_si128(tmp3, tmp2);
    tmp6 = _mm_xor_si128(tmp6, tmp3);
    return tmp6;
}


// Gueron and Kounavis fig. 8
__m128i gfmul_fig8(__m128i H1, __m128i H2, __m128i H3, __m128i H4, __m128i X1, __m128i X2, __m128i X3, __m128i X4 )
{
    __m128i H1_X1_lo, H1_X1_hi,H2_X2_lo, H2_X2_hi, H3_X3_lo, H3_X3_hi, H4_X4_lo, H4_X4_hi, lo, hi;
    __m128i tmp0, tmp1, tmp2, tmp3;
    __m128i tmp4, tmp5, tmp6, tmp7;
    __m128i tmp8, tmp9;
    H1_X1_lo = _mm_clmulepi64_si128(H1, X1,0x00);
    H2_X2_lo = _mm_clmulepi64_si128(H2, X2,0x00);
    H3_X3_lo = _mm_clmulepi64_si128(H3, X3, 0x00);
    H4_X4_lo = _mm_clmulepi64_si128(H4, X4,0x00);
    lo = _mm_xor_si128(H1_X1_lo, H2_X2_lo);
    lo = _mm_xor_si128(lo, H3_X3_lo);
    lo = _mm_xor_si128(lo, H4_X4_lo);
    H1_X1_hi = _mm_clmulepi64_si128(H1, X1, 0x11);
    H2_X2_hi = _mm_clmulepi64_si128(H2, X2, 0x11);
    H3_X3_hi = _mm_clmulepi64_si128(H3, X3, 0x11);
    H4_X4_hi = _mm_clmulepi64_si128(H4, X4,0x11);
    hi = _mm_xor_si128(H1_X1_hi, H2_X2_hi);
    hi = _mm_xor_si128(hi, H3_X3_hi);
    hi = _mm_xor_si128(hi, H4_X4_hi);
    tmp0 = _mm_shuffle_epi32(H1, 78);
    tmp4 = _mm_shuffle_epi32(X1, 78);
    tmp0 = _mm_xor_si128(tmp0, H1);
    tmp4 = _mm_xor_si128(tmp4, X1);
    tmp1 = _mm_shuffle_epi32(H2, 78);
    tmp5 = _mm_shuffle_epi32(X2, 78);
    tmp1 = _mm_xor_si128(tmp1, H2);
    tmp5 = _mm_xor_si128(tmp5, X2);
    tmp2 = _mm_shuffle_epi32(H3, 78);

    tmp6 = _mm_shuffle_epi32(X3, 78);
    tmp2 = _mm_xor_si128(tmp2, H3);
    tmp6 = _mm_xor_si128(tmp6, X3);
    tmp3 = _mm_shuffle_epi32(H4, 78);
    tmp7 = _mm_shuffle_epi32(X4, 78);
    tmp3 = _mm_xor_si128(tmp3, H4);
    tmp7 = _mm_xor_si128(tmp7, X4);
    tmp0 = _mm_clmulepi64_si128(tmp0,tmp4, 0x00);
    tmp1 = _mm_clmulepi64_si128(tmp1,tmp5, 0x00);
    tmp2 = _mm_clmulepi64_si128(tmp2,tmp6, 0x00);
    tmp3 = _mm_clmulepi64_si128(tmp3,tmp7, 0x00);
    tmp0 = _mm_xor_si128(tmp0, lo);
    tmp0 = _mm_xor_si128(tmp0, hi);
    tmp0 = _mm_xor_si128(tmp1, tmp0);
    tmp0 = _mm_xor_si128(tmp2, tmp0);
    tmp0 = _mm_xor_si128(tmp3, tmp0);
    tmp4 = _mm_slli_si128(tmp0, 8);
    tmp0 = _mm_srli_si128(tmp0, 8);
    lo = _mm_xor_si128(tmp4, lo);
    hi = _mm_xor_si128(tmp0, hi);
    tmp3 = lo;
    tmp6 = hi;
    tmp7 = _mm_srli_epi32(tmp3, 31);
    tmp8 = _mm_srli_epi32(tmp6, 31);
    tmp3 = _mm_slli_epi32(tmp3, 1);
    tmp6 = _mm_slli_epi32(tmp6, 1);
    tmp9 = _mm_srli_si128(tmp7, 12);
    tmp8 = _mm_slli_si128(tmp8, 4);
    tmp7 = _mm_slli_si128(tmp7, 4);
    tmp3 = _mm_or_si128(tmp3, tmp7);
    tmp6 = _mm_or_si128(tmp6, tmp8);
    tmp6 = _mm_or_si128(tmp6, tmp9);
    tmp7 = _mm_slli_epi32(tmp3, 31);
    tmp8 = _mm_slli_epi32(tmp3, 30);
    tmp9 = _mm_slli_epi32(tmp3, 25);
    tmp7 = _mm_xor_si128(tmp7, tmp8);
    tmp7 = _mm_xor_si128(tmp7, tmp9);
    tmp8 = _mm_srli_si128(tmp7, 4);
    tmp7 = _mm_slli_si128(tmp7, 12);
    tmp3 = _mm_xor_si128(tmp3, tmp7);
    tmp2 = _mm_srli_epi32(tmp3, 1);
    tmp4 = _mm_srli_epi32(tmp3, 2);
    tmp5 = _mm_srli_epi32(tmp3, 7);
    tmp2 = _mm_xor_si128(tmp2, tmp4);
    tmp2 = _mm_xor_si128(tmp2, tmp5);
    tmp2 = _mm_xor_si128(tmp2, tmp8);
    tmp3 = _mm_xor_si128(tmp3, tmp2);
    tmp6 = _mm_xor_si128(tmp6, tmp3);
    return tmp6;
}

// simple but fast implementation of GHASH (may not match exact spec. but should
// be just as fast)
__m128i GHASH_m128(__m128i key, const uint64_t * string64,
                   const size_t length) {
    unsigned int i = 0;
    __m128i answer = _mm_setzero_si128();
    size_t lengthm128 = length / 2;
    const __m128i * string = (const __m128i *) string64;

    if(lengthm128 > 4) {
        __m128i key2 = gfmul_fig5(key, key);
        __m128i key3 = gfmul_fig5(key2, key);
        __m128i key4 = gfmul_fig5(key2, key2);
        for(; i+4 < lengthm128; i+= 4) {
            __m128i B1 = _mm_loadu_si128(string + i );
            __m128i B2 = _mm_loadu_si128(string + i + 1);
            __m128i B3 = _mm_loadu_si128(string + i + 2);
            __m128i B4 = _mm_loadu_si128(string + i + 3);
            answer = _mm_xor_si128(answer,B1);
            answer = gfmul_fig8(key, key2, key3, key4,  B4, B3, B2, answer);
        }
    }
    for(; i < lengthm128; i++) {
        __m128i B = _mm_loadu_si128(string + i );
        answer =  gfmul_fig5(key,  _mm_xor_si128(answer,B));
    }
    if(length & 1) {// if odd
        __m128i B = _mm_loadl_epi64(string + lengthm128);
        answer =  gfmul_fig5(key,  _mm_xor_si128(answer,B));
    }
    return answer;
}


// for testing purposes, we need a 64-bit version of GHASH so we
// compute the full 128-bit version and then we just keep the least significant bits
// rs should point to 128 random bits.
uint64_t GHASH64bit(const void* rs, const uint64_t * string,
                    const size_t length) {
    __m128i key = _mm_loadu_si128((__m128i*) rs );
    __m128i answer = GHASH_m128(key, string, length);
    return _mm_cvtsi128_si64(answer);
}
#endif /* GHASH_H_ */
