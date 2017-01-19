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

static inline void mulBoth(uint64_t * x, uint64_t * y) {
  unsigned __int128 X = *x;
  unsigned __int128 Y = *y;
  X = X * Y;
  *x = X;
  *y = X >> 64;
}

static inline void addCarry(uint64_t * a, uint64_t * b, 
                            uint64_t c, uint64_t d,
                            uint64_t e, uint64_t f) {

  *a = c + e;
  __asm__ ("adcq %1, %2" : "+r"(d): "r"(f));
  *b = d;
  
  // __asm__ ("addq    %r8, %rdx \n\t"
  //          "adcq    %r9, %rcx \n\t"
  //          "movq    %rdx, (%rdi) \n\t"
  //          "movq    %rcx, (%rsi)");
           // : "=r"(*a), "=r"(*b)
           // : "r"(c), "r"(d), "r"(e), "r"(f));
}

// void addCarry(uint64_t * a, uint64_t * b, uint64_t c, uint64_t d) {
//   __asm__ ("addq    (%rdi), %rdx \n"
//         "adcq    (%rsi), %rcx \n"
//         "movq    %rdx, (%rdi) \n"
//            "movq    %rcx, (%rsi)");
  
// }

// static inline void mulBoth(uint64_t * x, uint64_t * y) {
//   uint64_t lo, hi;
//   __asm__("mulq %3" : "=a,a"(lo), "=d,d"(hi) : "%0,0"(x), "r,m"(y));
//   return hi;
// }

// Universal hashing of x and y. deltaDietz, like most hash functions
// in this file, relies on the following fact, referenced in the
// Badger paper as well as "A Construction Method for Optimally
// Universal Hash Families and its Consequences for the Existence of
// RBIBDs": if H is epsilon-almost delta universal, then the family
// consisting of the functions h'(x,y) = x + h(y), where h is from H,
// is universal.
//
// The proof is short: let (x0,y0) != (x1,y1). Then Pr_h(x0 + h(y0) =
// x1 + h(y1)) = Pr_h(x0 - x1 = h(y0) - h(y1)), by arithmetic. If y0 =
// y1, then x0 != x1 and the probability reduces to Pr_h(x0 = x1) =
// 0. Otherwise, by the delta-universality of H, the probability of
// h(y0) - h(y1) being ANY constant, INCLUDING x0 - x1, is low.
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
static inline __m256i clUniv512(const __m256i r128x2, const __m256i d0, const __m256i d1) {
  __m256i result = _mm256_xor_si256(r128x2, d0);

  // I think cast gets the lower bits of an __m256i
  __m128i x = _mm256_castsi256_si128(result);
  __m128i y = _mm256_extracti128_si256(result, 1);

  x = _mm_clmulepi64_si128(x, x, 1);
  y = _mm_clmulepi64_si128(y, y, 1);

  // I believe cast puts into the lower order bits of an __m256i
  result = _mm256_castsi128_si256(x);
  result = _mm256_inserti128_si256(result, y, 1);

  return _mm256_xor_si256(result, d1);
}

// In the interest of reducing code duplication, it is useful to
// parameterize treehashes by the primitive they use. By "primitive",
// I mean the building block that hashes R^2 down to R for some
// universe R. or instance, deltaDietz hashes uint64_t^2 down to
// uint64_t.
//
// We'll need a few ingredients:
struct MultiplyShift {
  // First, we need to know what R is in the notation above. This is
  // called an Atom, because from the point of view of the
  // parameterized hash algorithm, it cannot be broken up into
  // component parts.
  typedef uint64_t Atom;

 private:
  // Second, we need to know the type of the random data needed to
  // hash Atom^2 down to Atom.
  typedef ui128 Rand;

 public:
  // Third, we need to know the size of Atom in bytes. The sizeof
  // keyword won't always be enough is clients of this class, since
  // sometimes Atom is an array and "decays" into pointers.
  const static size_t ATOM_SIZE = sizeof(uint64_t);
  // Fourth, we need a way to copy Atoms.
  inline static void AtomCopy(Atom *x, const Atom &y) { *x = y; }
  // Fifth, we need to expose the hash family itself.
  inline void Hash(Atom *out, int i, const Atom &in0, const Atom &in1) const {
    *out = deltaDietz(r[i], in0, in1);
  }

  // Sixth, we need to store a particular array of random values. We
  // have to store these because some primitives need to preprocess
  // the random bits into an array of Rands. There is an example
  // below, but an easy one to think of is hashing modulo a prime. For
  // those hash families, we get data in, say, 8-bit chunks, but we
  // need random values modulo p, so we have to pre-process.
  //
  // We will need one Rand for each level of the tree hash. When we
  // pre-process, these values may require extra storage space. Since
  // that depends on the primitive, we keep the data in this struct
  // itself.
 private:
  const Rand *r;

 public:
  // In a primitive that needs preprocessing, depth indicates how many
  // random items to take - it is the depth of the treehash. We take
  // as an argument a const void ** that we update to point to random
  // data that this part of the tree hash will not need, so that it is
  // fresh when it is needed at the end of the tree hash.
  explicit MultiplyShift(const void **rvoid, const size_t depth)
      : r(reinterpret_cast<const Rand *>(*rvoid)) {
    *rvoid = reinterpret_cast<const void *>(r + depth);
  }
};

struct NH {
  static const bool alignmentRequired = false;
  typedef NH Unaligned;
  typedef ui128 Atom;

 private:
  typedef ui128 Rand;

 public:
  const static size_t ATOM_SIZE = sizeof(ui128);
  inline static void AtomCopy(Atom *x, const Atom &y) {
    (*x)[0] = y[0];
    (*x)[1] = y[1];
  }
  inline void Hash(Atom *out, const int i, const Atom &in0,
                   const Atom &in1) const {
    // Use a temporary Atom in case out == &in0 or out == &in1
    Atom tmp0 = {in0[0] + r[i][0], in0[1] + r[i][1]};
    mulBoth(&tmp0[0], &tmp0[1]);
    addCarry(&out[0][0], &out[0][1], tmp0[0], tmp0[1], in1[0], in1[1]);
  }

  explicit NH(const void **rvoid, const size_t depth)
      : r(reinterpret_cast<const Rand *>(*rvoid)) {
    *rvoid = reinterpret_cast<const void *>(r + depth);
  }

  inline static uint64_t Reduce(const void**rvoid, const Atom& x) {
    const ui128 * r128 = *reinterpret_cast<const ui128 **>(rvoid);
    *rvoid = reinterpret_cast<const void *>(r128+1);
    return deltaDietz(*r128, x[0], x[1]);
  }

 private:
  const Rand *r;
};

struct NHCLunaligned {
  static const bool alignmentRequired = false;
  typedef NHCLunaligned Unaligned;
  typedef __m128i Atom;

 private:
  typedef __m128i Rand;

 public:
  const static size_t ATOM_SIZE = sizeof(__m128i);
  inline static void AtomCopy(Atom *x, const Atom &y) { *x = y; }
  inline void Hash(Atom *out, const int i, const Atom &in0,
                   const Atom &in1) const {
    const Atom rk = r[i];
    const Atom tmp0 = _mm_lddqu_si128(&in0);
    const Atom tmp1 = _mm_clmulepi64_si128(rk, tmp0, 0x00);
    const Atom tmp2 = _mm_clmulepi64_si128(rk, tmp0, 0x11);
    const Atom tmp3 = _mm_xor_si128(tmp1, tmp2);
    *out = _mm_xor_si128(tmp3, _mm_lddqu_si128(&in1));
  }

  explicit NHCLunaligned(const void **rvoid, const size_t depth)
      : r(reinterpret_cast<const Rand *>(*rvoid)) {
    *rvoid = reinterpret_cast<const void *>(r + depth);
  }

  inline static uint64_t Reduce(const void **rvoid, const Atom &x) {
    const ui128 *r128 = *reinterpret_cast<const ui128 **>(rvoid);
    *rvoid = reinterpret_cast<const void *>(r128 + 1);
    uint64_t tmp[2];
    memcpy(tmp, &x, sizeof(x));
    return deltaDietz(*r128, tmp[0], tmp[1]);
  }

 private:
  const Rand *r;
};

// Basic mltilinear, with two multiplications and one addition:
struct NHCL {
  static const bool alignmentRequired = true;
  typedef NHCLunaligned Unaligned;
  typedef __m128i Atom;

 private:
  typedef __m128i Rand;

 public:
  const static size_t ATOM_SIZE = sizeof(__m128i);
  inline static void AtomCopy(Atom *x, const Atom &y) { *x = y; }
  inline void Hash(Atom *out, const int i, const Atom &in0,
                   const Atom &in1) const {
    const Atom rk = r[i];
    const Atom tmp1 = _mm_clmulepi64_si128(rk, in0, 0x00);
    const Atom tmp2 = _mm_clmulepi64_si128(rk, in0, 0x11);
    const Atom tmp3 = _mm_xor_si128(tmp1, tmp2);
    *out = _mm_xor_si128(tmp3, in1);
  }

  explicit NHCL(const void **rvoid, const size_t depth)
      : r(reinterpret_cast<const Rand *>(*rvoid)) {
    *rvoid = reinterpret_cast<const void *>(r + depth);
  }

  inline static uint64_t Reduce(const void **rvoid, const Atom &x) {
    const ui128 *r128 = *reinterpret_cast<const ui128 **>(rvoid);
    *rvoid = reinterpret_cast<const void *>(r128 + 1);
    uint64_t tmp[2];
    memcpy(tmp, &x, sizeof(x));
    return deltaDietz(*r128, tmp[0], tmp[1]);
  }

 private:
  const Rand *r;
};

struct CLNHunaligned {
  static const bool alignmentRequired = false;
  typedef CLNHunaligned Unaligned;
  typedef __m128i Atom;

 private:
  typedef __m128i Rand;

 public:
  const static size_t ATOM_SIZE = sizeof(__m128i);
  inline static void AtomCopy(Atom *x, const Atom &y) { *x = _mm_lddqu_si128(&y); }
  inline void Hash(Atom *out, const int i, const Atom &in0,
                   const Atom &in1) const {
    Atom tmp = _mm_xor_si128(_mm_lddqu_si128(&in0), r[i]);
    tmp = _mm_clmulepi64_si128(tmp, tmp, 1);
    tmp = _mm_xor_si128(tmp, _mm_lddqu_si128(&in1));
    *out = tmp;
  }

  explicit CLNHunaligned(const void **rvoid, const size_t depth)
      : r(reinterpret_cast<const Rand *>(*rvoid)) {
    *rvoid = reinterpret_cast<const void *>(r + depth);
  }

  inline static uint64_t Reduce(const void**rvoid, const Atom& x) {
    const ui128 * r128 = *reinterpret_cast<const ui128 **>(rvoid);
    *rvoid = reinterpret_cast<const void *>(r128+1);
    uint64_t tmp[2];
    memcpy(tmp, &x, sizeof(x));
    return deltaDietz(*r128, tmp[0], tmp[1]);
  }

 private:
  const Rand *r;
};


struct CLNH {
  static const bool alignmentRequired = true;
  typedef CLNHunaligned Unaligned;
  typedef __m128i Atom;

 private:
  typedef __m128i Rand;

 public:
  const static size_t ATOM_SIZE = sizeof(__m128i);
  inline static void AtomCopy(Atom *x, const Atom &y) { *x = y; }
  inline void Hash(Atom *out, const int i, const Atom &in0,
                   const Atom &in1) const {
    Atom tmp = _mm_xor_si128(in0, r[i]);
    tmp = _mm_clmulepi64_si128(tmp, tmp, 1);
    tmp = _mm_xor_si128(tmp, in1);
    *out = tmp;
  }

  explicit CLNH(const void **rvoid, const size_t depth)
      : r(reinterpret_cast<const Rand *>(*rvoid)) {
    *rvoid = reinterpret_cast<const void *>(r + depth);
  }

  inline static uint64_t Reduce(const void**rvoid, const Atom& x) {
    const ui128 * r128 = *reinterpret_cast<const ui128 **>(rvoid);
    *rvoid = reinterpret_cast<const void *>(r128+1);
    uint64_t tmp[2];
    memcpy(tmp, &x, sizeof(x));
    return deltaDietz(*r128, tmp[0], tmp[1]);
  }

 private:
  const Rand *r;
};

struct CLNHx2 {
  typedef __m256i Atom;

  const static size_t ATOM_SIZE = sizeof(__m256i);

  inline static void AtomCopy(Atom *x, const Atom &y) { *x = y; }

  explicit CLNHx2(const void **rvoid, const size_t depth) {
    const __m128i *const r128 = reinterpret_cast<const __m128i *>(*rvoid);
    for (size_t i = 0; i < depth; ++i) {
      r[i] = _mm256_broadcastsi128_si256(r128[i]);
    }
    *rvoid = reinterpret_cast<const void *>(r128 + depth);
  }

  inline void Hash(Atom *out, const int i, const Atom &in0,
                   const Atom &in1) const {
    *out = clUniv512(r[i], in0, in1);
  }

 private:
  typedef __m256i Rand;
  Rand r[62];
};

struct NHavx {
  static const bool alignmentRequired = true;
  typedef CLNHunaligned Unaligned;
  typedef __m256i Atom;

 private:
  typedef __m256i Rand;

 public:
  const static size_t ATOM_SIZE = sizeof(Atom);
  inline static void AtomCopy(Atom *x, const Atom &y) { *x = y; }
  inline void Hash(Atom *out, const int i, const Atom &in0,
                   const Atom &in1) const {
    Atom input = _mm256_add_epi32(in0, r[i]);
    const Atom hi = _mm256_srli_epi64(input, 32);
    input = _mm256_mul_epu32(input, hi);
    *out = _mm256_add_epi64(in1, input);
  }

  explicit NHavx(const void **rvoid, const size_t depth)
      : r(reinterpret_cast<const Rand *>(*rvoid)) {
    *rvoid = reinterpret_cast<const void *>(r + depth);
  }

  inline static uint64_t Reduce(const void **rvoid, const Atom &x) {
    int64_t y[4] = {_mm256_extract_epi64(x, 0), _mm256_extract_epi64(x, 1),
        _mm256_extract_epi64(x, 2), _mm256_extract_epi64(x, 3)};
    uint64_t *z = reinterpret_cast<uint64_t *>(y);
    const ui128 *r128 = *reinterpret_cast<const ui128 **>(rvoid);
    *rvoid = reinterpret_cast<const void *>(r128 + 2);
    z[0] = deltaDietz(*r128, z[0], z[1]);
    z[2] = deltaDietz(*r128, z[2], z[3]);
    return deltaDietz(*(r128 + 1), z[0], z[2]);
  }

 private:
  const Rand *r;
};

// This primitive just wraps another primitive, but operates on arrays
// of a given static size.
template <typename T, size_t n>
struct Wide {
  typedef typename T::Atom Atom[n];

  const static size_t ATOM_SIZE = T::ATOM_SIZE * n;

  inline static void AtomCopy(Atom *x, const Atom &y) {
    for (size_t i = 0; i < n; ++i) {
      T::AtomCopy(&(*x)[i], y[i]);
    }
  }

  explicit Wide(const void **rvoid, const size_t depth) : t(rvoid, depth) {}

  inline void Hash(Atom *out, const int j, const Atom &in0,
                   const Atom &in1) const {
    for (size_t i = 0; i < n; ++i) {
      t.Hash(&(*out)[i], j, in0[i], in1[i]);
    }
  }

 private:
  T t;
};

#endif
