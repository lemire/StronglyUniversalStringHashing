#ifndef SIMPLE_TREEHASH
#define SIMPLE_TREEHASH

#include <stdlib.h>
#include "util.hh"

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
  const uint64_t result = bigendian(*r128, readFrom[0], length);
  free(level);
  return result;
}

#endif
