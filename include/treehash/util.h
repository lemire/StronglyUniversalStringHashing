#ifndef TREEHASH_UTIL
#define TREEHASH_UTIL

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

#endif
