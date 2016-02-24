// CHeck that every word in an input string is used in the output, by
// varying that word and seeing that it changes the output.

#include <cassert>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string.h>

#include <iostream>

using namespace std;

extern "C" {
#include "hashfunctions64bits.h"
#include "pcg.h"
}

#include "treehash/binary-treehash.hh"
#include "treehash/generic-treehash.hh"

struct NamedFunc {
  const hashFunction64 f;
  const string name;
  NamedFunc(const hashFunction64& f, const string& name) : f(f), name(name) {}
};

#define NAMED(f) NamedFunc(f, #f)

NamedFunc hashFunctions[] = {
    NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, CLNH, 1>)),
    NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, NHCL, 7>)),
    NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, NH, 7>))};

const int HowManyFunctions64 =
    sizeof(hashFunctions) / sizeof(hashFunctions[0]);

void testunused() {
  printf("[%s] %s\n", __FILE__, __func__);
  assert(HowManyFunctions64 <= 64);
  int lengthStart = 1, lengthEnd = 2048; // inclusive
  int i;
  int length;
  uint64_t randbuffer[150] __attribute__ ((aligned (16)));// 150 should be plenty

  uint64_t * intstring;
  // We need 32 bytes of alignment for working with __m256i's
  if (posix_memalign((void **)(&intstring), 32, sizeof(uint64_t)*lengthEnd)) {
    cerr << "Failed to allocate " << lengthEnd << " words." << endl;
    return 1;
  }
  for (i = 0; i < 150; ++i) {
      randbuffer[i] = pcg64_random();
  }
  for (i = 0; i < lengthEnd; ++i) {
      intstring[i] = pcg64_random();
  }
  for (length = lengthStart; length <= lengthEnd; length += 1) {
    cout << length << " ";
    for (int place = 0; place < length; ++place) {
      for (i = 0; i < HowManyFunctions64; ++i) {
          const hashFunction64 thisfunc64 = hashFunctions[i].f;
          const auto first_run = thisfunc64(randbuffer, intstring, length);
          const auto old_val = intstring[place];
          const uint64_t new_val = old_val + 1; //rand() || ((uint64_t)(rand()) << 32);
          intstring[place] = new_val;
          const auto second_run = thisfunc64(randbuffer, intstring, length);
          if (first_run == second_run) {
            cerr << "uhoh: " << length << " " << place << " " << hashFunctions[i].name << endl
                 << first_run << " " << second_run << endl
                 << old_val << " " << new_val << endl;
            const auto old_val2 = intstring[place];
            const uint64_t new_val2 = old_val2 + 1; //rand() || ((uint64_t)(rand()) << 32);
            intstring[place] = new_val2;
            const auto second_run2 = thisfunc64(randbuffer, intstring, length);
            if (second_run2 == second_run) {
              cerr << second_run << " " << second_run2 << endl
                   << new_val << " " << new_val2 << endl;
            }
          }
      }
      }
      //printf("\n");
  }
  free(intstring);
}

int main(int c, char ** arg) {
  (void) (c);
  (void) (arg);

}
