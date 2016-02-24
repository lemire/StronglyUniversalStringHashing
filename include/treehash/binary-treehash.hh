#ifndef BINARY_TREEHASH
#define BINARY_TREEHASH

#include "util.hh"

// This version of treehash is based on the recursive treehash. That
// version was basically an depth-first traversal in which first the
// left subtree was traversed, then the right, then finally the root.
//
// The overhead of the function calls in the recursive treehash can be
// be reduced by maintaining the stack ourselves. This can be done by
// analogy to incrementing a binary number.
//
// As a reminder, here ais a tree used to calculate the hash of the string x0,...,x7
//
//  x0    x1  x2    x3  x4    x5  x6    x7
//    \  /      \  /      \  /      \  /
//     y0__    __y1        y2__    __y3
//         \  /                \  /
//          z0_______    _______z1
//                   \  /
//                    w0
//
// In recursive_treehash, the nodes are traversed in the order:
//
// x0,x1,y0,x2,x3,y1,z0,x4,x5,y2,x6,x7,y7,z1,w0
//
// Let's write these out as paths from the root, where L denotes taking
// a left branch and R a right one, and _ denotes the empty path:
//
// x0: LLL
// x1: LLR
// y0: LL
// x2: LRL
// x3: LRR
// y1: LR
// z0: L
// x4: RLL
// x5: RLR
// y2: RL
// x6: RRL
// x7: RRR
// y3: RR
// z1: R
// w0: _
//
// Now we'll write a series 3-bit binary numbers. They will represent
// the sequence of numbers 0..7,0, and they will show the process of
// doing the carries that happen when the number is
// incremented. During the carry process, the numbers are temporarily
// allowed one 2-bit.
//
// 0: 000
// 1: 001
// 2: 002 -> 010
// 3: 011
// 4: 012 -> 020 -> 100
// 5: 101
// 6: 102 -> 110
// 7: 111
// 0: 112 -> 120 -> 200 -> 000
//
// Now write these binary numbers as strings over {'L', 'R'}, with 0
// mapping to L and 1 to R, dropping any 2 digits and andy lower-order
// digits after 2s. This gives us the same list of tree paths as above.
//
// This connection between trees and binary numbers is well known. For
// instance, one can use it to show that insertions in a 2-3 tree need
// to rebalance modify O(1) nodes, amortized, and it is also used in
// "finger trees" where more exotic numeral systems can make inserts
// and deletes rebalance O(1) nodes in the worst case.
//
// We can use it to keep track of our location in a tree traversal
// using a small amount of data. First, we need to store all of the
// non-leaf values at the nodes we are actively working on. These are
// the partial hash values we are agglomerating into our final answer.
//
// Second, we need to maintain our path in the tree. This can be done
// with an unsigned int that doubles as our counter of how many xi's
// we have read. Note that we only need to store a partial hash value
// for levels where the path turns right, since those are the levels
// where the left value has already been calculated. This corresponds
// to 1's in the binary representation of the path.
//
// Third, we need to maintain the carry digit. The location can be
// inferred from the int counter representing the path, but we also
// need to store a partial hash value representing the carry digit
// being propagated.


// Perform the carry operation, moving `carry` up the path
// (prepresented by `workspace`. This function performs exactly
// `limit` carries.
static inline void cascade(const ui128* r128, uint64_t workspace[64],
                           const uint64_t carry, int limit) {
  const uint64_t* last = &carry;
  int i = 0;
  for (; i+1 < limit; ++i) {
    workspace[i] = deltaDietz(r128[i], workspace[i], *last);
    last = &workspace[i];
  }
  workspace[i+1] = deltaDietz(r128[i], workspace[i], *last);
}

// Once the partial calculations are complete, we need to combine them
// in binary trees that are non-complete. This function does so. It
// uses randomness only if count is not a power of 2. The return vale
// is the location in `workspace` of the collected value.
static inline size_t rollup(const ui128* r128, uint64_t workspace[64],
                            size_t count) {
  // We need a place to start, so we look for the least significant
  // 1-bit in count.
  int i = __builtin_ctzll(count);
  size_t last = i;
  ++i;
  count >>= i;
  // While there are any 1-bits left, we combine the previous 1-bit with
  // the next one.
  for (; count > 0; ++i, count >>= 1) {
    if (count & 1) {
      workspace[i] = deltaDietz(r128[i], workspace[i], workspace[last]);
      last = i;
    }
  }
  return last;
}

uint64_t binary_treehash(const void* rvoid, const uint64_t* data,
                         const size_t length) {
  const ui128* r128 = (const ui128*)rvoid;
  uint64_t workspace[64];
  for (size_t i = 0; i < length; ++i) {
    if (i & 1) {
      cascade(r128, workspace, data[i], __builtin_ctzll(i + 1));
    } else {
      workspace[0] = data[i];
    }
  }
  const size_t last = rollup(r128, workspace, length);
  r128 = &r128[(length & (length-1)) ? (last+1) : last];
  return bigendian(*r128, workspace[last], length);
}

// This is just linke binary_treehash, but works on generic primitives
// like CLNH. See util.hh, which has examples that show the interface
// that T must provide. Perhaps most important are Atom and Hash. Atom
// is the type that T operates on, and Hash reduces two Atoms to one,
// almost-universally.
//
// Not threadsafe.
template <typename T>
struct GenericBinaryTreehash {
 private:

  // We store a value of type T (rather than just using static
  // functions) because it may need to construct its randomness. See
  // util.hh.
  T primitive;
  typedef typename T::Atom Atom;

 public:

  // As in the primitives in util.hh, the constructor adjusts *r so it
  // does not point at any randomness to be used by this object.
  GenericBinaryTreehash(const void** r, const int depth)
      : primitive(r, depth){};


  // The method treehash is not static because it needs to use
  // primitive and workspace.
  Atom* treehash(const Atom* data, const size_t length) {
    for (size_t i = 0; i < length; ++i) {
      if (i & 1) {
        cascade(data[i], __builtin_ctzll(i + 1));
      } else {
        T::AtomCopy(&workspace[0], data[i]);
      }
    }
    const size_t last = rollup(length);
    return &workspace[last];
  }

 private:
  Atom workspace[64];

  inline void cascade(const Atom& carry, int limit) {
    const Atom* last = &carry;
    int i = 0;
    for (; i + 1 < limit; ++i) {
      primitive.Hash(&workspace[i], i, workspace[i], *last);
      last = &workspace[i];
    }
    primitive.Hash(&workspace[i + 1], i, workspace[i], *last);
  }

  inline size_t rollup(size_t count) {
    int i = __builtin_ctzll(count);
    size_t last = i;
    ++i;
    count >>= i;
    for (; count > 0; ++i, count >>= 1) {
      if (count & 1) {
        primitive.Hash(&workspace[i], i, workspace[i], workspace[last]);
        last = i;
      }
    }
    return last;
  }

};

// GenericBinaryTreeHash uses T::AtomCopy one out of every two
// iterations of the main loop of treehash. In fact, it copies the
// entire input. ZeroCopyGenericBinaryTreehash performs the same
// calculation without the copies.
//
// Not threadsafe.
template <typename T>
struct ZeroCopyGenericBinaryTreehash {
 private:
  T primitive;

  typedef typename T::Atom Atom;

  // ZeroCopyGenericBinaryTreehash has some member variables that
  // GenericBinaryTreeHash does not.
  //
  // First, it stores the depth that was passed in in the
  // constructor. This will be used later to simplify rollup.
  int depth;
  // Second, it uses overflow as a storage location to put the result
  // of T::Hash when workspace[0] is full.
  Atom overflow;
  // Additionally, workspace can be a bit smaller, since the levels
  // are all going to be bumped up. In GenericBinaryTreeHash,
  // workspace[0] was the input, workspace[1] was the t::Hash of the
  // input, workspace[2] was the T::Hash of workspace[1], and so
  // in. In this class, the input is not stored in a member
  // variable. overflow and workspace[0] are hashes of the input, and
  // so on.
  Atom workspace[63];

  inline void cascade(int limit) {
    const Atom* last = &overflow;
    int i = 0;
    for (; i + 1 < limit; ++i) {
      primitive.Hash(&workspace[i], i+1, workspace[i], *last);
      last = &workspace[i];
    }
    primitive.Hash(&workspace[i + 1], i+1, workspace[i], *last);
  }

  // The parameter extra represents any data from the input that didn'
  // have a pair, so it used only in rollup, not cascade or the
  // treehash main loop. It is repurposed in the function body to
  // represent the Atom most recently produced.
  inline void rollup(const Atom * extra, size_t count) {
    int i = __builtin_ctzll(count);
    size_t last = i;
    ++i;
    count >>= i;
    if (!extra) extra = &workspace[last];
    for (; count > 0; ++i, count >>= 1) {
      if (count & 1) {
        primitive.Hash(&workspace[i], i+1, workspace[i], *extra);
        last = i;
        extra = &workspace[last];
      }
    }
  }

 public:
  Atom* treehash(const Atom* data, const size_t length) {
    for (size_t i = 0; i+2 <= length; i += 2) {
      if (i & 2) {
        // workspace[0] is full
        primitive.Hash(&overflow, 0, data[i], data[i+1]);
        cascade(__builtin_ctzll(i + 1));
      } else {
        primitive.Hash(&workspace[0], 0, data[i], data[i+1]);
      }
    }
    rollup((length & 1) ? &data[length-1] : nullptr, length/2);
    // Since length >= 2, depth >= 1, so workspace[depth-1] has been populated.
    return &workspace[depth-1];
  }

  ZeroCopyGenericBinaryTreehash(const void** r, const int depth)
    : primitive(r, depth), depth(depth) {};
};

// ZeroCopyGenericBinaryTreehash iterated over the input two Atoms at
// once. BoostedZeroCopyGenericBinaryTreehash increases that two four
// by reusing workspace[0] and workspace[1] as scratch space.
//
// Not threadsafe.
template <typename T>
struct BoostedZeroCopyGenericBinaryTreehash {
 private:
  T primitive;

  typedef typename T::Atom Atom;

  int depth;
  Atom overflow;
  Atom workspace[63];

  inline void cascade(int limit) {
    const Atom* last = &overflow;
    int i = 0;
    for (; i + 1 < limit; ++i) {
      primitive.Hash(&workspace[i], i+1, workspace[i], *last);
      last = &workspace[i];
    }
    primitive.Hash(&workspace[i + 1], i+1, workspace[i], *last);
  }

  inline size_t rollup(const Atom * extra, size_t count) {
    if (extra) {
      int i = __builtin_ctzll(count);
      primitive.Hash(&workspace[i], i+1, workspace[i], *extra);
      int last = i;
      ++i;
      count >>= i;
      for (; count > 0; ++i, count >>= 1) {
        if (count & 1) {
          primitive.Hash(&workspace[i], i+1, workspace[i], workspace[last]);
          last = i;
        }
      }
      return last;
    } else {
      int i = __builtin_ctzll(count);
      int last = i;
      ++i;
      count >>= i;
      for (; count > 0; ++i, count >>= 1) {
        if (count & 1) {
          primitive.Hash(&workspace[i], i+1, workspace[i], workspace[last]);
          last = i;
        }
      }
      return last;
    }
  }

 public:
  Atom* treehash(const Atom* data, const size_t length) {
    size_t i = 0;
    for (;i+8 <= length; i += 8) {
      // workspace[0] and workspace[1] are empty.
      primitive.Hash(&workspace[0], 0, data[i], data[i+1]);
      primitive.Hash(&workspace[1], 0, data[i+2], data[i+3]);
      primitive.Hash(&workspace[1], 1, workspace[0], workspace[1]);
      // workspace[1] is full, but overflow and workspace[0] are empty.
      primitive.Hash(&workspace[0], 0, data[i+4], data[i+5]);
      primitive.Hash(&overflow, 0, data[i+6], data[i+7]);
      cascade(__builtin_ctzll((i + 8)/2));
    }
    if (i+4 <= length) {
      if (i & 4) {
        // workspace[1] is full, but overflow and workspace[0] are empty.
        primitive.Hash(&workspace[0], 0, data[i], data[i+1]);
        primitive.Hash(&overflow, 0, data[i+2], data[i+3]);
        cascade(__builtin_ctzll((i + 4)/2));
      } else {
        // workspace[0] and workspace[1] are empty.
        primitive.Hash(&workspace[0], 0, data[i], data[i+1]);
        primitive.Hash(&workspace[1], 0, data[i+2], data[i+3]);
        primitive.Hash(&workspace[1], 1, workspace[0], workspace[1]);
      }
      i += 4;
    }
    if (i+2 <= length) {
      if (i & 2) {
        primitive.Hash(&overflow, 0, data[i], data[i+1]);
        cascade(__builtin_ctzll(i + 1));
      } else {
        primitive.Hash(&workspace[0], 0, data[i], data[i+1]);
      }
    }
    const size_t last = rollup((length & 1) ? &data[length-1] : nullptr, length/2);
    return &workspace[last];
  }

  BoostedZeroCopyGenericBinaryTreehash(const void** r, const int depth)
    : primitive(r, depth), depth(depth) {};
};

#endif
