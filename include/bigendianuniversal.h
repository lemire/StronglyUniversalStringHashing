#ifndef BIGENDIANUNIVERSAL_H
#define BIGENDIANUNIVERSAL_H

// A family of hash functions from U to D is said to be
// "epsilon-almost big-endian universal" when, for any x != y in U and
// s < lg D, the probability that h(x) >> s = h(y) >> s is less than
// 2^s * epsilon. These functions may be suitable for use in hash
// tables when the most significant bits, rather than the least
// significant bits, are used as the index into the table.
//
// The famous multiply-shift family from Z_{2^u} to Z_{2^d} of
// { h_a : a \in Z_{2^u}, a odd } where h_a(x) is defined as (a * x) >> (u-d)
// is 2^{1-d}-almost big-endian universal.
//
// By iterating this hash function over a string with words in
// Z_{2^u}, we get an L*2^{1-d}-almost big-endian universal family,
// where L is the maximum length of a string.

#include <string.h>

// Unsigned 128-bit integers
typedef struct U128 {
  uint64_t hi, lo;
} u128;

// Multuply two uint64_ts and get the most significant 64 bits
// of the result; Intel only.
static inline uint64_t hi64mul(uint64_t x, uint64_t y) {
  uint64_t lo, hi;
  __asm__ ("mulq %3" : "=a,a" (lo), "=d,d" (hi) : "%0,0" (x), "r,m" (y));
  return hi;
}

// Multiply two u128s, but don't compute the least significant 64 bits
// or the most significant 128 bits.
static inline u128 multHi128(u128 x, u128 y) {
  const u128 result = {x.hi * y.lo + x.lo * y.hi + hi64mul(x.lo, y.lo),
                       0};
  return result;
}

// An L*2^{1-d} almost bigendian universal hash function on strings of
// 64-bit words, producing a 64-bit output.
uint64_t hornerHash(const void * randomSource,
                    const uint64_t * x,
                    const size_t length) {
  u128 h;
  memcpy(&h, randomSource, sizeof(u128));
  // h must be odd:
  h.lo |= 1;
  // We treat the length as the first word in the string, ensuring
  // that no string is a prefix of any other.
  u128 accum = {length, x[0]};
  accum = multHi128(h, accum);
  for (size_t i = 1; i < length; ++i) {
    // accum.hi holds the hash value we have accumulated so far. We
    // put the next word to hash into accum.lo to make one two-word
    // integer that we hash with h:
    accum.lo = x[i];
    accum = multHi128(h, accum);
  }
  return accum.hi;
}

// Hashes two 64-bit words (newer and accum) down to one,
// universally. h0 and h1 must be chosen uniformly at random.
//
// This method first rehashes `accum` using `h0` and `h1`. It is
// essentially (accum * (h0 + (h1<<64))) >> 64, where all operations
// are performed on 128-bit words. This is a variant on the
// multiply-shift hashing in Dietzfelbinger's "Universal hashing and
// k-wise independent random variables via integer arithmetic without
// primes". It doesn't use the addition component. This makes it only
// delta-universal, not strongly universal. However, that is enough so
// that adding `newer` maintains the universailty, as has been noted
// by many.
static inline void univHash(const uint64_t h0, const uint64_t h1,
                            const uint64_t newer, uint64_t *accum) {
  *accum = newer + *accum * h1 + hi64mul(*accum, h0);
}

// Horner's method can only dispatch one 128-bit multiplication at a
// time, since each loop iteration depends on the one before
// it. unrolledHorner changes the order in which the words are hashed
// and can hash multiple words simultaneously, if the processor
// supports it.
#define DECLARE_UNROLLED_HORNER(MANY)                                        \
  uint64_t unrolledHorner##MANY(const void *randomSource, const uint64_t *x, \
                                const size_t length) {                     \
    uint64_t accums[MANY];                                                   \
    accums[0] = length;                                                      \
    accums[1] = x[1];                                                        \
    for (size_t i = 2; i < MANY; ++i) {                                      \
      accums[i] = (i < length) ? x[i] : 0;                                   \
    }                                                                        \
    size_t i = (length < MANY) ? length : MANY;                              \
    const uint64_t *r64 = (const uint64_t *)randomSource;                    \
    for (; i + MANY <= length; i += MANY) {                                  \
      for (size_t j = 0; j < MANY; ++j) {                                    \
        univHash(r64[0], r64[1], x[i + j], &accums[j]);                      \
      }                                                                      \
    }                                                                        \
    for (size_t j = 0; i+j < length; j += 1) {                               \
      univHash(r64[0], r64[1], x[i+j], &accums[j]);                          \
    }                                                                        \
    for (size_t j = 1; j < MANY; ++j) {                                      \
      univHash(r64[0], r64[1], accums[j], &accums[0]);                       \
    }                                                                        \
    return accums[0] * (r64[2] | ((uint64_t)1));                             \
  }

DECLARE_UNROLLED_HORNER(3)
DECLARE_UNROLLED_HORNER(4)
DECLARE_UNROLLED_HORNER(5)
DECLARE_UNROLLED_HORNER(6)
DECLARE_UNROLLED_HORNER(7)
DECLARE_UNROLLED_HORNER(8)
DECLARE_UNROLLED_HORNER(9)

// One other way to calculate a 64-bit hash value is to calculate two
// 32-bit hash values. In order to make this almost big-endian
// universal, we have to rehash these two 32-bit values with an
// almost-big-endian universal function, too.
//
// This function is further from being big-endian universal, in that
// its epsilon is larger. The two 32-bit hash values each have
// collision probability approximately L*2^{-31}, so the probability
// that they both collide is about L^2 * 2^{-62}, rather than
// L*{2^-63} in the constructions above. This method is thus better
// suited for shorter strings than longer ones.
//
// This function is unrolled in a manenr similar to unrolledHorner,
// above.
uint64_t twiceHorner32(const void * randomSource,
                       const uint64_t * x,
                       const size_t length) {
  u128 h;
  memcpy(&h, randomSource, sizeof(u128));
  // h.lo and h.hi must be odd:
  h.lo |= 1;
  h.hi |= 1;
  if (1 == length) {
    const u128 tmp = {length, x[0]};
    return multHi128(h, tmp).hi;
  }
  if (2 == length) {
    u128 tmp = {length, x[0]};
    tmp = multHi128(h, tmp);
    tmp.lo = x[1];
    return multHi128(h, tmp).hi;
  }
  u128 accums[4];
  accums[0].hi = length; accums[1].hi = x[0];
  accums[2].hi = x[1]; accums[3].hi = x[2];
  size_t i = 3;
  // This is the main loop.
  for (; i + 3 < length; i += 4) {
    for (size_t j = 0; j < 4; ++j) {
      accums[j].lo = x[i+j] * h.lo;
      accums[j].hi *= h.hi;
      accums[j].hi &= 0xffffffff00000000ull;
      accums[j].hi |= accums[j].lo >> 32;
    }
  }
  // We might have 1, 2, or 3 words left over at the end that we
  // couldn't handle in our unrolled loop which could only do 4 at
  // once:
  for(; i < length; i += 1) {
    for (size_t j = 0; j < 1; ++j) {
      accums[j].lo = x[i+j] * h.lo;
      accums[j].hi *= h.hi;
      accums[j].hi &= 0xffffffff00000000ull;
      accums[j].hi |= accums[j].lo >> 32;
    }
  }
  // Finally, we combine all of the hash values we have already computed.
  for (size_t j = 1; j < 4; ++j) {
    accums[0].lo = accums[j].hi;
    accums[0] = multHi128(h, accums[0]);
  }
  return accums[0].hi;
}

// Fast universal hashing using the CLNH family:

// Using the random words in `r64`, universally hash `data` down to 64 bits
static inline uint64_t bigHashDown(const uint64_t *r64, const __m128i *data) {
  int64_t d64[2] = {_mm_extract_epi64(*data, 0),
                    _mm_extract_epi64(*data, 1)};
  uint64_t * u64 = (uint64_t *)d64;
  univHash(r64[0], r64[1], u64[0], &u64[1]);
  return u64[1];
}

// Big-endian universally hash `data` with the random bits in `r64`.
static inline uint64_t bie(const uint64_t r64, const uint64_t data) {
  return data * (r64 | 0x1ul);
}

// Hash universally `accum` and `data` using the random bits in `r128`.
static inline void clCombine(const __m128i r128, __m128i *accum,
                             const __m128i data) {
  *accum = _mm_xor_si128(*accum, r128);
  *accum = _mm_clmulepi64_si128(*accum, *accum, 1);
  *accum = _mm_xor_si128(*accum, data);
}

// Helper functions for hashing an array of __m128i `accum` and a
// single __m128i `extra` down into `accum[0]` using the random bits
// in `r128`.
static inline void clArrayCombineExtra2(const __m128i r128, __m128i *accum,
                                        const __m128i extra) {
  clCombine(r128, &accum[0], extra);
  clCombine(r128, &accum[0], accum[1]);
}

static inline void clArrayCombineExtra3(const __m128i r128, __m128i *accum,
                                        const __m128i extra) {
  clCombine(r128, &accum[0], extra);
  clCombine(r128, &accum[1], accum[2]);
  clCombine(r128, &accum[0], accum[1]);
}

static inline void clArrayCombineExtra4(const __m128i r128, __m128i *accum,
                                        const __m128i extra) {
  for (size_t i = 0; i < 2; ++i) {
    clCombine(r128, &accum[i], accum[i + 2]);
  }
  clArrayCombineExtra2(r128, accum, extra);
}

static inline void clArrayCombineExtra5(const __m128i r128, __m128i *accum,
                                        const __m128i extra) {
  clCombine(r128, &accum[0], extra);
  for (size_t i = 1; i <= 2; ++i) {
    clCombine(r128, &accum[i], accum[i + 2]);
  }
  clArrayCombineExtra2(r128, accum, accum[2]);
}

static inline void clArrayCombineExtra6(const __m128i r128, __m128i *accum,
                                        const __m128i extra) {
  for (size_t i = 0; i < 3; ++i) {
    clCombine(r128, &accum[i], accum[i + 3]);
  }
  clArrayCombineExtra3(r128, accum, extra);
}

static inline void clArrayCombineExtra8(const __m128i r128, __m128i *accum,
                                        const __m128i extra) {
  for (size_t i = 0; i < 4; ++i) {
    clCombine(r128, &accum[i], accum[i + 4]);
  }
  clArrayCombineExtra4(r128, accum, extra);
}

static inline void clArrayCombineExtra9(const __m128i r128, __m128i *accum,
                                        const __m128i extra) {
  clCombine(r128, &accum[0], extra);
  for (size_t i = 1; i <= 4; ++i) {
    clCombine(r128, &accum[i], accum[i + 4]);
  }
  clArrayCombineExtra4(r128, accum, accum[4]);
}

static inline void clArrayCombineExtra10(const __m128i r128, __m128i *accum,
                                         const __m128i extra) {
  for (size_t i = 0; i < 5; ++i) {
    clCombine(r128, &accum[i], accum[i + 5]);
  }
  clArrayCombineExtra5(r128, accum, extra);
}

static inline void clArrayCombineExtra11(const __m128i r128, __m128i *accum,
                                         const __m128i extra) {
  clCombine(r128, &accum[0], extra);
  for (size_t i = 1; i <= 5; ++i) {
    clCombine(r128, &accum[i], accum[i + 5]);
  }
  clArrayCombineExtra5(r128, accum, accum[5]);
}

static inline void clArrayCombineExtra12(const __m128i r128, __m128i *accum,
                                         const __m128i extra) {
  for (size_t i = 0; i < 6; ++i) {
    clCombine(r128, &accum[i], accum[i + 6]);
  }
  clArrayCombineExtra6(r128, accum, extra);
}

// Iterated hashing using the CLNH family
#define ITERATECL(MANY)                                                     \
  uint64_t iterateCL##MANY(const void *randomSource, const uint64_t *x,     \
                           const size_t length) {                           \
    const uint64_t *r64 = &(((const uint64_t *)(randomSource))[2]);         \
    const __m128i *z = (const __m128i *)x;                                  \
    const size_t zlen = length / 2;                                         \
    __m128i accum[MANY];                                                    \
    for (size_t j = 0; j < MANY; ++j) {                                     \
      accum[j] = (j < zlen) ? _mm_lddqu_si128(&z[j]) : _mm_setzero_si128(); \
    }                                                                       \
    size_t i = (zlen >= MANY) ? MANY : zlen;                                \
    const __m128i r128 = _mm_lddqu_si128((const __m128i *)randomSource);    \
    for (; i + MANY <= zlen; i += MANY) {                                   \
      for (size_t j = 0; j < MANY; ++j) {                                   \
        clCombine(r128, &accum[j], _mm_lddqu_si128(&z[i + j]));             \
      }                                                                     \
    }                                                                       \
    for (size_t j = 0; j < MANY; ++j) {                                     \
      if (i + j < zlen) {                                                   \
        clCombine(r128, &accum[j], _mm_lddqu_si128(&z[i + j]));             \
      }                                                                     \
    }                                                                       \
    clArrayCombineExtra##MANY(                                              \
        r128, accum,                                                        \
        _mm_set_epi64x(length, (length & 1) ? x[length - 1] : 0));          \
    const uint64_t bhd = bigHashDown(r64, &accum[0]);                       \
    return bie(r64[2], bhd);                                                \
  }

ITERATECL(8)
ITERATECL(9)
ITERATECL(10)
ITERATECL(11)
ITERATECL(12)

#undef ITERATECL

static inline __m128i clCombineFar(const __m128i r128, const __m128i x,
                                   const __m128i y) {
  __m128i result = _mm_xor_si128(x, r128);
  result = _mm_clmulepi64_si128(result, result, 1);
  result = _mm_xor_si128(result, y);
  return result;
}

// Unrolled Carter & Wegman tree-hash with clCombineFar as the
// reducing function. This is basically "Badger - A Fast and Provably
// Secure MAC", by Boesgaard et al.

#define HALVE(MANY)                                                       \
  static inline void halve##MANY(const __m128i r128, const __m128i *from, \
                                 __m128i *to) {                           \
    for (size_t i = 0; i < MANY; ++i) {                                   \
      to[i] = clCombineFar(r128, from[2 * i], from[2 * i + 1]);           \
    }                                                                     \
  }

#define HALVE_LOAD(MANY)                                                      \
  static inline void halveLoad##MANY(const __m128i r128, const __m128i *from, \
                                     __m128i *to) {                           \
    for (size_t i = 0; i < MANY; ++i) {                                       \
      to[i] = clCombineFar(r128, _mm_lddqu_si128(&from[2 * i]),               \
                           _mm_lddqu_si128(&from[2 * i + 1]));                \
    }                                                                         \
  }

#define TREECL(MANY)                                                        \
  uint64_t treeCL##MANY(const void *randomSource, const uint64_t *x,        \
                        const size_t length) {                              \
    const __m128i *r128 = (const __m128i *)randomSource;                    \
    const size_t depth = 64 - __builtin_clzll(1 + length);                  \
    __m128i rLevel[64];                                                     \
    for (size_t i = 0; i < depth; ++i) {                                    \
      rLevel[i] = _mm_lddqu_si128(&r128[i]);                                \
    }                                                                       \
    __m128i tree[64][2 * MANY];                                             \
    size_t fill[64];                                                        \
    for (size_t i = 0; i < depth; ++i) {                                    \
      fill[i] = 0;                                                          \
    }                                                                       \
    const __m128i *z = (const __m128i *)x;                                  \
    const size_t zlen = length / 2;                                         \
    size_t i = 0;                                                           \
    for (; i + 2 * MANY <= zlen; i += 2 * MANY) {                           \
      for (size_t j = 0; 2 * MANY == fill[j]; ++j) {                        \
        halve##MANY(rLevel[j + 1], tree[j], &tree[j + 1][fill[j + 1]]);     \
        fill[j] = 0;                                                        \
        fill[j + 1] += MANY;                                                \
      }                                                                     \
      halveLoad##MANY(rLevel[0], &z[i], &tree[0][fill[0]]);                 \
      fill[0] += MANY;                                                      \
    }                                                                       \
    size_t max_fill_level = depth - 1;                                      \
    for (; fill[max_fill_level] > 0; --max_fill_level) {                    \
    }                                                                       \
    for (size_t j = 0; 2 * MANY == fill[j]; ++j) {                          \
      halve##MANY(rLevel[j + 1], tree[j], &tree[j + 1][fill[j + 1]]);       \
      fill[j] = 0;                                                          \
      fill[j + 1] += MANY;                                                  \
    }                                                                       \
    for (; i < zlen; i += 2) {                                              \
      tree[0][fill[0]] = clCombineFar(rLevel[0], _mm_lddqu_si128(&z[i]),    \
                                      _mm_lddqu_si128(&z[i + 1]));          \
      ++fill[0];                                                            \
    }                                                                       \
    const __m128i final =                                                   \
        _mm_set_epi64x(length, (length & 1) ? x[length - 1] : 0);           \
    tree[0][fill[0]] =                                                      \
        (i < zlen) ? clCombineFar(rLevel[0], _mm_lddqu_si128(&z[i]), final) \
                   : final;                                                 \
    ++fill[0];                                                              \
    i = 0;                                                                  \
    for (; (i < max_fill_level) || (fill[i] > 1); ++i) {                    \
      size_t j = 0;                                                         \
      for (; j + 2 <= fill[i]; j += 2) {                                    \
        tree[i + 1][fill[i + 1]] =                                          \
            clCombineFar(rLevel[i + 1], tree[i][j], tree[i][j + 1]);        \
        ++fill[i + 1];                                                      \
      }                                                                     \
      if (j < fill[i]) {                                                    \
        tree[i + 1][fill[i + 1]] = tree[i][j];                              \
        ++fill[i + 1];                                                      \
      }                                                                     \
    }                                                                       \
    const uint64_t *r64 = (const uint64_t *)randomSource;                   \
    r64 += 2 * MANY;                                                        \
    const uint64_t bhd = bigHashDown(r64, &tree[i][0]);                     \
    return bie(r64[2], bhd);                                                \
  }

#define TREECLALL(MANY) HALVE(MANY) HALVE_LOAD(MANY) TREECL(MANY)

TREECLALL(8)
TREECLALL(9)
TREECLALL(10)

#endif  // BIGENDIANUNIVERSAL_H
