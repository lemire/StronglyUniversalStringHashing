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

uint64_t simple_treehash(const void * rvoid, const uint64_t * data, const size_t length) {
  const ui128* r128 = (const ui128 *)rvoid;
  uint64_t * level = (uint64_t *)malloc(sizeof(uint64_t) * (length+1)/2);
  const uint64_t * readFrom = data;
  for (size_t lengthLeft = length; lengthLeft > 1; lengthLeft = (lengthLeft+1)/2) {
    size_t i = 0;
    for (; i+1 < lengthLeft; i += 2) {
      level[i/2] = deltaDietz(*r128, readFrom[i], readFrom[i+1]);
    }
    if (lengthLeft & 1) {
      level[i/2] = readFrom[i];
    }
    readFrom = level;
    ++r128;
  }
  const uint64_t result = bigendian(*r128, level[0], length);
  free(level);
  return result;
}
