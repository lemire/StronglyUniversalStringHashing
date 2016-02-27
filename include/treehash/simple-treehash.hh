// This file shows simple variants on tree hashing, following Carter
// and Wegman's "Universal Classes of Hash Functions" and Boesgaard et
// al.'s "Badger - A Fast and Provably Secure MAC".

#ifndef SIMPLE_TREEHASH
#define SIMPLE_TREEHASH

#include <stdlib.h>
#include "util.hh"

// Hash a string universally, but don't hash in its length. This is
// only almost universal for domains containing strings of only one
// length.
//
// The first argument, r128, it incremented so that at return time
// *r128 points to random bits that have not yet been used.
static inline uint64_t simple_treehash_without_length(const ui128 **r128,
                                                      const uint64_t *data,
                                                      const size_t length,
                                                      uint64_t *const level) {
  const uint64_t *readFrom = data;
  for (size_t lengthLeft = length; lengthLeft > 1;
       lengthLeft = (lengthLeft + 1) / 2) {
    size_t i = 0;
    for (; i + 1 < lengthLeft; i += 2) {
      level[i / 2] = deltaDietz(**r128, readFrom[i], readFrom[i + 1]);
    }
    if (lengthLeft & 1) {
      level[i / 2] = readFrom[i];
    }
    readFrom = level;
    ++*r128;
  }
  return readFrom[0];
}

// This is the same idea as simple_treehash_without_length, but
// genericized and with a reduced number of copies. It also does not
// update r.
template <typename T>
static inline const typename T::Atom *generic_simple_treehash_without_length(
    const T &prefill, const typename T::Atom *data, const size_t length,
    typename T::Atom *level) {
  // length >= 2, level has size (length + 1)/2
  if ((length & 1) and (data != level)) {
    T::AtomCopy(&level[0], data[0]);
  }
  const typename T::Atom *readFrom = data;
  int level_idx = 0;
  for (size_t lengthLeft = length; lengthLeft > 1;
       lengthLeft = (lengthLeft + 1) / 2) {
    size_t i = lengthLeft & 1;
    for (; i + 2 <= lengthLeft; i += 2) {
      prefill.Hash(&level[(i + 1) / 2], level_idx, readFrom[i],
                   readFrom[i + 1]);
    }
    readFrom = level;
    ++level_idx;
  }
  return &readFrom[0];
}

// The same as simple_treehash_without_length, but first consolidates
// two different pieces of input data (`lhs` and `data`) into one.
template <typename T>
static inline uint64_t split_simple_treehash_without_length(
    const ui128 **r128, const typename T::Atom &lhs, const uint64_t *data,
    const size_t length) {
  uint64_t level[2 * T::ATOM_SIZE];
  memcpy(level, &lhs, T::ATOM_SIZE);
  memcpy(reinterpret_cast<char *>(level) + T::ATOM_SIZE, data,
         length * sizeof(uint64_t));
  const uint64_t result =
    simple_treehash_without_length(r128, level, T::ATOM_SIZE + length, level);
  return result;
}

template <typename T, size_t N>
static inline const typename T::Atom *
split_generic_simple_treehash_without_length(const void **rvoid,
                                             typename T::Atom level[],
                                             const typename T::Atom lhs[],
                                             const uint64_t *data,
                                             const size_t length) {
  typedef typename T::Atom Atom;
  memcpy(level, lhs, N * T::ATOM_SIZE);
  memcpy(reinterpret_cast<char *>(level) + N * T::ATOM_SIZE, data,
         length * sizeof(uint64_t));
  memset(reinterpret_cast<char *>(level) + N * T::ATOM_SIZE +
             length * sizeof(uint64_t),
         0, N * T::ATOM_SIZE - length * sizeof(uint64_t));
  T prefill(rvoid, 64 - __builtin_clzll(N));
  return generic_simple_treehash_without_length(
      prefill, level,
      N + (sizeof(uint64_t) * length + sizeof(Atom) - 1) / sizeof(Atom), level);
}

uint64_t simple_treehash(const void *rvoid, const uint64_t *data,
                         const size_t length) {
  const ui128 *r128 = (const ui128 *)rvoid;
  uint64_t *const level =
      reinterpret_cast<uint64_t *>(malloc(sizeof(uint64_t) * (length + 1) / 2));
  const uint64_t result =
      simple_treehash_without_length(&r128, data, length, level);
  free(level);
  return bigendian(*r128, result, length);
}

// This is the same as simple_treehash, but uses a stack-allocated
// workspace for performance reasons.
template <size_t n>
uint64_t short_simple_treehash(const void *rvoid, const uint64_t *data,
                               const size_t length) {
  const ui128 *r128 = (const ui128 *)rvoid;
  uint64_t level[(n + 1) / 2];
  const uint64_t result =
      simple_treehash_without_length(&r128, data, length, level);
  return bigendian(*r128, result, length);
}

// This is like simple_treehash, but works on generic hashing
// primitives (see util.hh).
template <typename T>
uint64_t generic_simple_treehash(const void *rvoid, const uint64_t *data,
                                 const size_t length) {
  typedef typename T::Atom Atom;
  static const size_t ATOM_WORD_SIZE = T::ATOM_SIZE / sizeof(uint64_t);
  static_assert(sizeof(uint64_t) * ATOM_WORD_SIZE == T::ATOM_SIZE,
                "sizeof(Atom) must be a multiple of sizeof(uint64_t)");
  if (length < 2 * ATOM_WORD_SIZE) {
    return short_simple_treehash<2 * ATOM_WORD_SIZE>(rvoid, data, length);
  }
  const size_t atom_length = length / ATOM_WORD_SIZE;
  const Atom *atom_data = reinterpret_cast<const Atom *>(data);
  // We will reduce atom_length Atoms down to 1 Atom. This requires
  // ceiling(log2(atom_length)) levels of tree hashing. For x > 1,
  // __builtin_clzll(x) is 64 - ceiling(log2(x)), unless x is a power
  // of 2, in which case it is 65 - ceiling(log2(x)). subtracting 1
  // from x only changes its __builtin_clzll if x is a power of 2. In
  // that case, it increases __builtin_clzll by 1. Thus, 64 -
  // __builtin_clzll(length-1) is ceiling(log2(length)).
  const size_t levels_count = 64 - __builtin_clzll(atom_length - 1);

  // // We need levels_count Rands; rnext points to one past the end of
  // // that. We will use it when the tree hashing part of this function
  // // is complete.
  // const void *rnext = reinterpret_cast<const void *>(
  //     reinterpret_cast<const typename T::Rand *>(rvoid) + levels_count);

  // Load the randomness into a hashing object
  T prefill(&rvoid, levels_count);
  // We need workspace that is 32-byte aligned to work with types like
  // __m256i. AVX512 will need 64-byte aligned data, presumably.
  Atom *level;
  if (posix_memalign((void **)(&level), 32,
                     T::ATOM_SIZE * ((atom_length + 1) / 2))) {
    exit(1);
  }
  const Atom *tree_result = generic_simple_treehash_without_length<T>(
      prefill, atom_data, atom_length, level);
  const size_t data_read = ATOM_WORD_SIZE * atom_length;
  const ui128 *r128 = reinterpret_cast<const ui128 *>(rvoid);
  const uint64_t result = split_simple_treehash_without_length<T>(
      &r128, *tree_result, &data[data_read], length - data_read);
  free(level);
  return bigendian(*r128, result, length);
}

// In simple_cl_treehash below, we will need to use 128 bits of
// randomness duplicated into one __m256i.
static inline void prefill_rand128x2(__m256i * r128x2, const __m128i * r128,
                                     const size_t n) {
  for (size_t i = 0; i < n; ++i) {
    r128x2[i] = _mm256_broadcastsi128_si256(r128[i]);
  }
}

uint64_t simple_cl_treehash(const void * rvoid, const uint64_t * data,
                            const size_t length) {
  if (length < 8) return short_simple_treehash<8>(rvoid, data, length);
  const __m128i* r128 = (const __m128i *)rvoid;
  // We will reduce length * 64 bits down to 256 bits. This requires
  // ceiling(log2(length)) - 2 units of randomness, where one unit is
  // 128 bits. For x > 1, __builtin_clzll(x) is 64 - ceiling(log2(x)),
  // unless x is a power of 2, in which case it is 65 -
  // ceiling(log2(x)). subtracting 1 from x only changes its
  // __builtin_clzll if x is a power of 2. In that case, it increases
  // __builtin_clzll by 1. Thus, 64 - __builtin_clzll(length-1) is
  // ceiling(log2(length)).
  const size_t cl_levels = 62 - __builtin_clzll(length-1);
  __m256i r128x2[cl_levels];
  prefill_rand128x2(r128x2, r128, cl_levels);
  // r128x2cursor tracks where the next randomness to use is. Once we
  // have used up r128x2, we will start using r128[cl_levels] (and
  // later indices as needed).
  __m256i * r128x2cursor = r128x2;
  // Moving from (data, length) in uint64_t-land to (d256, length256)
  // in __m256i-land.
  const __m256i * d256 = reinterpret_cast<const __m256i *>(data);
  const size_t length256 = length/4;
  // We need workspace that is 32-byte aligned to work with __m256i
  __m256i * level;
  if (posix_memalign((void **)(&level), 32,
                     sizeof(__m256i) * (length256 + 1) / 2)) {
    exit(1);
  }
  const __m256i *readFrom = d256;
  for (size_t lengthLeft = length256; lengthLeft > 1;
       lengthLeft = (lengthLeft+1)/2) {
    size_t i = 0;
    for (; i+1 < lengthLeft; i += 2) {
      level[i / 2] = clUniv512(*r128x2cursor, readFrom[i], readFrom[i + 1]);
    }
    if (lengthLeft & 1) {
      level[i/2] = readFrom[i];
    }
    readFrom = level;
    ++r128x2cursor;
  }
  // We've hashed most of the string down to a sinlge __m256i. What's left:
  //
  // 1. const __m256i level[0]
  //
  // 2. const uint64_t data[4*length256] .. data[length-1], at most 3 of them
  //
  // 3. const size_t length
  //
  // 4. const r128[levels] .. r128[63]
  //
  // We put #1 and #2 into a single buffer and handle them using
  // simple_treehash.
  uint64_t final_level[7];
  _mm256_storeu_si256(reinterpret_cast<__m256i *>(final_level), level[0]);
  free(level);
  size_t i = 0;
  for (;4*length256 + i< length; ++i) {
    final_level[4+i] = data[4*length256 + i];
  }
  const ui128 * ir128 = reinterpret_cast<const ui128*>(r128 + cl_levels);
  const uint64_t result =
      simple_treehash_without_length(&ir128, final_level, 4 + i, final_level);
  return bigendian(*ir128, result, length);
}

#endif
