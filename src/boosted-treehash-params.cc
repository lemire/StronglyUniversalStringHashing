// This program benchmarks different depths of boosting in
// boosted_treehash and prints tuples of
//
// 1. Length, in 64-bit words, of the strings being hashed
// 2. Percentile
// 3. The fastest parameter for that percentile
//
// For instance, the line "1023 74 8" means that when hashing strings
// of length 1023, the 72nd percentile boosted_treehash<8> isfaster
// than the 72nd percentile boosted_treehash<n> for other values of
// n that were tested;

#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

#include "treehash/boosted-treehash.hh"
extern "C" {
#include "timers.h"
}

void fill_random(uint64_t *data, size_t n) {
  const int rfd = open("/dev/urandom", O_RDONLY);
  if (-1 == rfd) {
    const int err = errno;
    fprintf(stderr, "%s\n", strerror(err));
    exit(err);
  }
  char *const cdata = (char *)data;
  for (size_t i = 0; i < n; i += read(rfd, &cdata[i], sizeof(uint64_t) * n - i))
    ;
  (void)close(rfd);
}

uint64_t *alloc_random(size_t n) {
  uint64_t *ans = (uint64_t *)malloc(sizeof(uint64_t) * n);
  if (!ans) {
    fprintf(stderr, "Failed allocation of %zd words\n", n);
    exit(1);
  }
  fill_random(ans, n);
  return ans;
}

// A dummy value we update only to force computation
uint64_t sum = 0;

// One more than the maximum parameter N to boosted_treehash<N>
const size_t MAX_DEPTH = 9;

template <size_t N>
void bench(vector<ticks> *out, const void *r, const uint64_t *data,
           const size_t i) {
  const ticks start = startRDTSC();
  sum += boosted_treehash<N>(r, data, i);
  const ticks finish = startRDTSC();
  out->push_back(finish - start);
  bench<N + 1>(out + 1, r, data, i);
}

template <>
void bench<MAX_DEPTH>(vector<ticks> *, const void *, const uint64_t *,
                      const size_t) {
  return;
}

int main() {
  size_t max_len = 2048;
  const uint64_t *const data = alloc_random(max_len);
  const void *const r64 = alloc_random(128);
  size_t iters = 10000;
  vector<vector<ticks> > cycles(MAX_DEPTH, vector<ticks>());
  size_t samples = 100;
  cout << "# length percentile optimal-treeboost " << endl;
  for (size_t i = 1; i < max_len; i = 1 + 1.001 * i) {
    for (size_t j = 0; j < MAX_DEPTH; ++j) {
      cycles[j].clear();
    }
    for (size_t j = 0; j < iters; ++j) {
      bench<0>(&cycles[0], r64, data, i);
    }
    for (size_t j = 0; j < MAX_DEPTH; ++j) {
      sort(cycles[j].begin(), cycles[j].end(), std::greater<ticks>());
    }
    for (size_t j = 0; j <= samples; ++j) {
      size_t loc = (iters * j) / samples;
      if (loc >= iters) loc = iters - 1;
      ticks min_val = numeric_limits<ticks>::max();
      int min_idx = -1;
      for (size_t k = 0; k < MAX_DEPTH; ++k) {
        if (cycles[k][loc] < min_val) {
          min_val = cycles[k][loc];
          min_idx = k;
        }
      }
      cout << i << " " << j << " " << min_idx << endl;
    }
    cout << endl;
  }
  if (0 == sum) {
    cerr << "# Magic happens: sum of all hashes is 0!" << endl;
    return 1;
  }
};
