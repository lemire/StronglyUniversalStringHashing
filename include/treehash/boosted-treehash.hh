#ifndef BOOSTED_TREEHASH
#define BOOSTED_TREEHASH

#include "util.hh"

// binary_treehash is still slower than simple_treehash. To avoid the
// bookkeeping of the carries from `cascade`, we can take the lowest k
// levels of the tree (equivalent: the full subtrees at level k) and
// hash them without testing any conditions - we can set at compile
// time which carries to perform. Additionally, if k is set at compile
// time, we do not need a loop and the jumps and conditionals in that.
//
// This essentially "boosts" the binary treehash up k levels. It is
// similar to the technique used in VHASH and CLHASH. In those
// functions, the lower level hash function (NH and CLNH,
// respectively) is fast, but requires a lot of randomness, while the
// higher level function (polynomial hashing) is slower but requires
// less randomness. Here. the lower level function is faster by using
// more temporary storage space and by using longer instruction space
// (the unrolled loop).

// We will still need `cascade` and `rollup` from binary treehash.
#include "binary-treehash.hh"
// To be able to set k at compile times for different values of k, we
// use template specialization.

// finish_treehash is a treehash that works like the simple treehash,
// but works on a mutable workspace and so does not need to allocate
// new space.
template <size_t log>
static inline uint64_t finish_treehash(const ui128* r128,
                                       uint64_t simple_workspace[1 << log]) {
  static const size_t power = 1ull << log;
  for (size_t i = 0; i < power / 2; ++i) {
    simple_workspace[i] =
        deltaDietz(*r128, simple_workspace[2 * i], simple_workspace[2 * i + 1]);
  }
  return finish_treehash<log - 1>(r128 + 1, simple_workspace);
}

template <>
inline uint64_t finish_treehash<0>(const ui128*, uint64_t simple_workspace[1]) {
  return simple_workspace[0];
}

// static_treehash is like finish_treehash, but cannot overwrite its
// input.
template <size_t log>
static inline uint64_t static_treehash(
    const ui128* r128, const uint64_t* data,
    uint64_t *simple_workspace) {
  static const size_t power = 1ull << log;
  for (size_t i = 0; i < power / 2; ++i) {
    simple_workspace[i] = deltaDietz(*r128, data[2 * i], data[2 * i + 1]);
  }
  return finish_treehash<log - 1>(r128 + 1, simple_workspace);
}

template <>
inline uint64_t static_treehash<0>(const ui128*, const uint64_t* data,
                                   uint64_t *) {
  return data[0];
}

// finalize_treehash<L> completes a treehash by hashing the words left
// over at the end of a string. The total length left over must be
// less than 1<<L. It puts some 1-bits in the lowest L slots of binary_workspace
template <size_t log>
static inline void finalize_treehash(
    const ui128* r128, const uint64_t* data, const size_t length,
    uint64_t binary_workspace[64], uint64_t * simple_workspace) {
  static const size_t power = 1ull << log;

  if (power / 2 <= length) {
    binary_workspace[log - 1] =
        static_treehash<log - 1>(r128, data, simple_workspace);
    finalize_treehash<log - 1>(r128, &data[power / 2], length - power / 2,
                               binary_workspace, simple_workspace);
  } else {
    finalize_treehash<log - 1>(r128, data, length, binary_workspace,
                               simple_workspace);
  }
}

template <>
inline void finalize_treehash<0>(const ui128*, const uint64_t*, const size_t,
                                 uint64_t*, uint64_t*) {}

template <size_t log>
uint64_t boosted_treehash(const void* rvoid, const uint64_t* data,
                          const size_t length) {
  const ui128* r128 = reinterpret_cast<const ui128*>(rvoid);
  uint64_t binary_workspace[64];
  static const size_t power = 1ull << log;
  uint64_t simple_workspace[power / 2];
  size_t i = 0;
  for (; i + power <= length; i += power) {
    if (i & power) {
      cascade(r128 + log, binary_workspace + log,
              static_treehash<log>(r128, &data[i], simple_workspace),
              __builtin_ctzll(i + power) - log);
    } else {
      binary_workspace[log] =
          static_treehash<log>(r128, &data[i], simple_workspace);
    }
  }
  finalize_treehash<log>(r128, &data[i], length - i, binary_workspace,
                         simple_workspace);
  const size_t last = rollup(r128, binary_workspace, length);
  r128 = &r128[(length & (length - 1)) ? (last + 1) : last];
  return bigendian(*r128, binary_workspace[last], length);
}

#endif
