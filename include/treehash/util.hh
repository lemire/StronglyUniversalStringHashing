#ifndef TREEHASH_UTIL
#define TREEHASH_UTIL

#include <cstdint>
#include "immintrin.h"

typedef uint64_t ui128[2];

// Multiply two unsigned 64-bit ints, producing a 128-bit int but
// returning only the high-order bits.
static inline uint64_t mulHi(uint64_t x, uint64_t y) {
  uint64_t lo, hi;
  __asm__("mulq %3" : "=a,a"(lo), "=d,d"(hi) : "%0,0"(x), "r,m"(y));
  return hi;
}
// Universal hashing of x and y
static inline uint64_t deltaDietz(const ui128 h, const uint64_t x,
                                   const uint64_t y) {
  return x + y * h[1] + mulHi(y, h[0]);
}

static inline uint64_t bigendian(const ui128 h, const uint64_t x,
                                 const uint64_t y) {
  const uint64_t hlo = h[0] | 0x1;
  return h[1] * x + hlo * y + mulHi(hlo, x);
}

// Universal hashing 512 bits down to 256 bits using 128 bits of
// random data, with collision probability 2^-64.
//
// r128x2 is 128 bits of random data, duplicated.
//
// The method we use is one that relies on:
//
// 1. CLNH
//
// 2. The fact that if H is a almost delta-universal family, we can
// get an almost universal family by using one value from H and then
// performing an addition in the codomain. This has been noted in the
// Badger paper, as well as by Woelfel in "A construction method for
// optimally universal hash families and its consequences for the
// existence of RBIBDs.".
static inline __m256i clUniv512(const __m256i r128x2, const __m256i d[2]) {
  __m256i result = _mm256_xor_si256(r128x2, d[0]);

  // I think cast gets the lower bits of an __m256i
  __m128i x = _mm256_castsi256_si128(result);
  __m128i y = _mm256_extracti128_si256(result, 1);

  x = _mm_clmulepi64_si128(x, x, 1);
  y = _mm_clmulepi64_si128(y, y, 1);

  // I believe cast puts into the lower order bits of an __m256i
  result = _mm256_castsi128_si256(x);
  result = _mm256_inserti128_si256(result, y, 1);

  return _mm256_xor_si256(result, d[1]);
}

#endif
