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

void fill_random(uint64_t * data, size_t n) {
  const int rfd = open("/dev/urandom", O_RDONLY);
  if (-1 == rfd) {
    const int err = errno;
    fprintf(stderr, "%s\n", strerror(err));
    exit(err);
  }
  char * const cdata = (char *)data;
  for (size_t i = 0; i < n; i += read(rfd, &cdata[i], sizeof(uint64_t)*n-i));
  (void)close(rfd);
}

uint64_t * alloc_random(size_t n) {
  uint64_t * ans = (uint64_t *)malloc(sizeof(uint64_t) * n);
  if (!ans) {
    fprintf(stderr, "Failed allocation of %zd words\n", n);
    exit(1);
  }
  fill_random(ans, n);
  return ans;
}

uint64_t sum = 0;

template<size_t N>
ticks bench(const void * r, const uint64_t * data, const size_t i) {
  const ticks start = startRDTSC();
  sum += boosted_treehash<N>(r, data, i);
  const ticks finish = startRDTSC();
  return finish-start;
}

int main() {
  size_t max_len = 2048;
  const uint64_t * const data = alloc_random(max_len);
  const void * const r64 = alloc_random(128);
  size_t iters = 10000;
  const size_t max_depth = 9;
  vector<vector<ticks> > cycles(max_depth, vector<ticks>(iters));
  size_t samples = 100;
  for (size_t i = 1; i < max_len; i = 1 + 1.01*i) {
    for (size_t j = 0; j < iters; ++j) {
      cycles[0][j] = bench<1>(r64, data, i);
      cycles[1][j] = bench<2>(r64, data, i);
      cycles[2][j] = bench<3>(r64, data, i);
      cycles[3][j] = bench<4>(r64, data, i);
      cycles[4][j] = bench<5>(r64, data, i);
      cycles[5][j] = bench<6>(r64, data, i);
      cycles[6][j] = bench<7>(r64, data, i);
      cycles[7][j] = bench<8>(r64, data, i);
      cycles[8][j] = bench<9>(r64, data, i);
    }
    for (size_t j = 0; j < max_depth; ++j) {
      sort(cycles[j].begin(), cycles[j].end(), std::greater<ticks>());
    }
    for (size_t j = 0; j <= samples; ++j) {
      size_t loc = (iters * j)/samples;
      if (loc >= iters) loc = iters-1;
      ticks min_val = numeric_limits<ticks>::max();
      int min_idx = -1;
      for (size_t k = 0; k < max_depth; ++k) {
        if (cycles[k][loc] < min_val) {
          min_val = cycles[k][loc];
          min_idx = k;
        }
      }
      /*
      for (size_t k = 0; k < max_depth; ++k) {
        min_val = min(min_val, cycles[k][loc]);
      }
      */
      cout << i << " " << j << " " << min_idx << endl;
    }
    cout << endl;
  }
  cerr << "# ignore " << sum << endl;
};
