#ifndef RECURSIVE_TREEHASH
#define RECURSIVE_TREEHASH

#include "util.hh"

// The basic treehash requires the allocation of additional space
// approximately half of the size of the input. This version avoids
// that, and allocates space proportional to the logarithm (base 2) of
// the input.
//
// The way it works is as follows: consider hashing a string of length
// 8. First, we allocate scratch space of size 4. We then halve the
// length of the string each time. Visually:
//
//  x0    x1  x2    x3  x4    x5  x6    x7
//    \  /      \  /      \  /      \  /
//     y0__    __y1        y2__    __y3
//         \  /                \  /
//          z0_______    _______z1
//                   \  /
//                    w0
//
// where yi is a hash of x(2*i) and x(2*i+1), using the first 128 bits
// in the random seed, and zi is a hash of y(2*i) and y(2*i+1) using
// the second 128 bits in the random seed, and so on.
//
// In simple_treehash, these values are calculated in the order
//
// y0,y1,y2,y3, z0,z1, w0
//
// Notice that after we calculate a value, we no longer use the values
// of its children in the tree. If we calculate the values in a
// different order, we can ensure that we need at most 1 + lg2(n)
// values (other than the input, x) at any time. Specifically, we can
// calculate in a greedy fashion -- after we calculate y0 and y1, we
// can calculate z0 and forget y0 and y1.
//
// To illustrate this, here are some diagrams of treess representing
// the calculation. The diagrams will represent a time series (or
// crude animation) of the calculation. In these diagrams, the values
// are no longer need to inspect will be elided, while the most recent
// values inspected or produces will have capital letters.
//
// First, for the simple treehash case:
//
//////////////////////////////////////////////////
//  X0    X1
//    \  /
//     Y0
//////////////////////////////////////////////////
//        X2    X3
//          \  /
//     y0    Y1
//////////////////////////////////////////////////
//              X4    X5
//                \  /
//     y0    y1    Y2
//////////////////////////////////////////////////
//                    X6    X7
//                      \  /
//     y0    y1    y2    Y3
//////////////////////////////////////////////////
//     Y0    Y1    y2    y3
//       \  /
//        Z0
//////////////////////////////////////////////////
//                 Y2    Y3
//                   \  /
//        z0          Z1
//////////////////////////////////////////////////
//          Z0_    _Z1
//             \  /
//              W0
//////////////////////////////////////////////////
//
// Note that as y3 is being calculated, y0, y1, and y2 are still
// resident in memory, as we will need hem to calculate the z
// level. In this case, that is only three extra cases, but, for a
// string of length 2*N, that is N y's that stay in memory.
//
// The greedy method looks like this:
//
//////////////////////////////////////////////////
//  X0    X1
//    \  /
//     Y0
//////////////////////////////////////////////////
//        X2    X3
//          \  /
//     y0    Y1
//////////////////////////////////////////////////
// (Note that these first two steps are the same as the simple
// method. Now that we have both y0 and y1, though, we can calculate z0)
//////////////////////////////////////////////////
//     Y0    Y1
//       \  /
//        Z0
//////////////////////////////////////////////////
//           X4    X5
//             \  /
//              Y2
//
//        z0
//////////////////////////////////////////////////
//                 X6    X7
//                   \  /
//              y2    Y3
//
//        z0
//////////////////////////////////////////////////
//             Y2    Y3
//               \  /
//        z0      Z1
//////////////////////////////////////////////////
//         Z0    Z1
//           \  /
//            W0
//////////////////////////////////////////////////
//
// This method needs at most a constant number of storage locations at
// each level of the tree. Space-wise, the simple treehash method is
// like breadth-first traversal while this greedy method is like
// post-order depth-first traversal.

// Hash a string (without hashing in the length) of length exactly
// (1 << depth). Uses `depth` random values from the array `r128`.
uint64_t perfect_treehash(const ui128* r128, const uint64_t* data,
                          const size_t depth) {
  if (0 == depth) return data[0];
  const uint64_t lhs = perfect_treehash(r128, data, depth - 1);
  const uint64_t rhs =
    perfect_treehash(r128, data + (1 << (depth - 1)), depth - 1);
  return deltaDietz(r128[depth - 1], lhs, rhs);
}

// Hash a string without hashing in the length.
uint64_t imperfect_treehash(const ui128* r128, const uint64_t* data,
                            const size_t length) {
  // depth is floor(log_2(length)).
  const size_t depth = 63 - __builtin_clzll(length);
  // This uses `depth` random values from r128.
  const uint64_t lhs = perfect_treehash(r128, data, depth);
  const size_t lhs_length = 1ull << depth;
  // If we can return early, we've used only `depth` random values
  // from `r128`. Otherwise, well have to use one more.
  if (length == lhs_length) return lhs;
  const uint64_t rhs =
      imperfect_treehash(r128, data + lhs_length, length - lhs_length);
  return deltaDietz(r128[depth], lhs, rhs);
}

uint64_t recursive_treehash(const void* rvoid, const uint64_t* data,
                            const size_t length) {
  const ui128* r128 = (const ui128*)rvoid;
  if (1 == length) return bigendian(r128[0], data[0], length);
  const uint64_t result = imperfect_treehash(r128, data, length);
  // This detects how many random items from `r128` were used to
  // calculate `result`.
  const size_t depth = 64 - __builtin_clzll(length-1);
  return bigendian(r128[depth], result, length);
}

#endif
