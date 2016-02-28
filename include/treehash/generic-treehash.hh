#ifndef GENERIC_TREEHASH
#define GENERIC_TREEHASH

#include "simple-treehash.hh"

// This is like simple_generic_treehash, but works on generic tree
// hashing primitives. Here, ALGO's constructor must use a logarithmic
// amount of randomness and ALGO::treehash must return a pointer to
// its result.
template <template<typename> class ALGO, typename T, size_t N>
uint64_t generic_treehash(const void *rvoid, const uint64_t *data,
                                 const size_t length) {
  if (T::alignmentRequired && (0 != (reinterpret_cast<size_t>(data) & (T::ATOM_SIZE - 1)))) {
    return generic_treehash<ALGO, typename T::Unaligned, N>(rvoid, data, length);
  }
  typedef typename Wide<T,N>::Atom Atom;
  static const size_t ATOM_WORD_SIZE = Wide<T,N>::ATOM_SIZE / sizeof(uint64_t);
  static_assert(sizeof(uint64_t) * ATOM_WORD_SIZE == Wide<T,N>::ATOM_SIZE,
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

  // The ALGO constructor moves rvoid forward past levels_count levels
  // of the randomness that T needs.
  ALGO<Wide<T,N> > hasher(&rvoid, levels_count);

  const Atom *tree_result = hasher.treehash(atom_data, atom_length);

  const size_t data_read = ATOM_WORD_SIZE * atom_length;
  typename T::Atom level[2*N];
  const typename T::Atom *result1 =
      split_generic_simple_treehash_without_length<T, N>(
          &rvoid, level, *tree_result, &data[data_read], length - data_read);
  const uint64_t result2 = T::Reduce(&rvoid, *result1);
  return bigendian(*reinterpret_cast<const ui128 *>(rvoid), result2, length);
}

#endif  // GENERIC_TREEHASH
