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
inline uint64_t hi64mul(uint64_t x, uint64_t y) {
  uint64_t lo, hi;
  __asm__ ("mulq %3" : "=a,a" (lo), "=d,d" (hi) : "%0,0" (x), "r,m" (y));
  return hi;
}

// Multiply two u128s, but don't compute the least significant 64 bits
// or the most significant 128 bits.
inline u128 multHi128(u128 x, u128 y) {
  return (u128) {.hi = x.hi * y.lo + x.lo * y.hi + hi64mul(x.lo, y.lo),
                 .lo = 0};
}

// An L*2^{1-d} almost strongly universal hash functions on strings of
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
  u128 accum = multHi128(h, (u128) {.hi = length, .lo = x[0]});
  for (size_t i = 1; i < length; ++i) {
    // accum.hi holds the hash value we have accumulated so far. We
    // put the next word to hash into accum.lo to make one two-word
    // integer that we hash with h:
    accum.lo = x[i];
    accum = multHi128(h, accum);
  }
  return accum.hi;
}

// Horner's method can only dispatch one 128-bit multiplication at a
// time, since each loop iteration depends on the one before
// it. unrolledHorner changes the order in which the words are hashed
// and can hash four words simultaneously, if the processor supports
// it.
uint64_t unrolledHorner(const void * randomSource,
                        const uint64_t * x,
                        uint64_t length) {
  u128 h;
  memcpy(&h, randomSource, sizeof(u128));
  // h must be odd:
  h.lo |= 1;
  if (1 == length) {
    return multHi128(h, (u128) {.hi = length, .lo = x[0]}).hi;
  }
  if (2 == length) {
    u128 tmp = multHi128(h, (u128) {.hi = length, .lo = x[0]});
    tmp.lo = x[1];
    return multHi128(h, tmp).hi;
  }
  u128 accums[4] = {(u128) {.hi = length}, (u128) {.hi = x[0]},
                    (u128) {.hi = x[1]}, (u128) {.hi = x[2]}};
  size_t i = 3;
  // This is the main loop:
  for (; i + 3 < length; i += 4) {
    for (size_t j = 0; j < 4; ++j) {
      accums[j].lo = x[i+j];
      accums[j] = multHi128(h, accums[j]);
    }
  }
  // We might have 1, 2, or 3 words left over at the end that we
  // couldn't handle in our unrolled loop which could only do 4 at
  // once:
  for(; i < length; i += 1) {
    for (size_t j = 0; j < 1; ++j) {
      accums[j].lo = x[i+j];
      accums[j] = multHi128(h, accums[j]);
    }
  }
  // Finally, we combine all of the hash values we have already computed.
  for (size_t j = 1; j < 4; ++j) {
    accums[0].lo = accums[j].hi;
    accums[0] = multHi128(h, accums[0]);
  }
  return accums[0].hi;
}

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
                        uint64_t length) {
  u128 h;
  memcpy(&h, randomSource, sizeof(u128));
  // h.lo and h.hi must be odd:
  h.lo |= 1;
  h.hi |= 1;
  if (1 == length) {
    return multHi128(h, (u128) {.hi = length, .lo = x[0]}).hi;
  }
  if (2 == length) {
    u128 tmp = multHi128(h, (u128) {.hi = length, .lo = x[0]});
    tmp.lo = x[1];
    return multHi128(h, tmp).hi;
  }
  u128 accums[4] = {(u128) {.hi = length}, (u128) {.hi = x[0]},
                    (u128) {.hi = x[1]}, (u128) {.hi = x[2]}};
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


#endif  // BIGENDIANUNIVERSAL_H
