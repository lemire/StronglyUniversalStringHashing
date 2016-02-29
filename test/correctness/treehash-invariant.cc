#include <cstdint>
#include <cstdio>
#include <cinttypes>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include "treehash/simple-treehash.hh"
#include "treehash/recursive-treehash.hh"
#include "treehash/binary-treehash.hh"
#include "treehash/boosted-treehash.hh"

bool fill_random(uint64_t * data, ssize_t n) {
  const int rfd = open("/dev/urandom", O_RDONLY);
  if (-1 == rfd) {
    fprintf(stderr, "%s\n", strerror(errno));
    return false;
  }
  char * const cdata = (char *)data;
  for (ssize_t i = 0; i < n; i += read(rfd, &cdata[i], sizeof(uint64_t)*n-i));
  return true;
}

int main(int argc, char ** argv) {
  ssize_t length = 1 << 15;
  if (argc > 1) {
    if (argc > 2) {
      fprintf(stderr, "Only one argument is premitted: the length\n");
      return 1;
    }
    if (sscanf(argv[1], "%zd", &length) != 1) {
      fprintf(stderr, "Argument \"%s\" is not an integer\n", argv[1]);
      return 1;
    }
  }
  uint64_t r[128];
  if (!fill_random(r, sizeof(r)/sizeof(r[0]))) {
    return 1;
  }
  uint64_t * data = nullptr;
  for (ssize_t i = 1; i <= length; ++i) {
    if (!(i & (i-1))) {
      printf("Success up to %zd\n", i-1);
      if (data) free(data);
      const ssize_t data_len = 2*i;
      data = reinterpret_cast<uint64_t *>(malloc(sizeof(uint64_t) * data_len));
      if (!data) {
        fprintf(stderr, "Malloc failed: %zd\n", data_len);
        return 1;
      }
      if (!fill_random(data, data_len)) return 1;
    }
    const uint64_t x[4] = {simple_treehash(r, data, i),
                           recursive_treehash(r, data, i),
                           binary_treehash(r, data, i),
                           boosted_treehash<5>(r, data, i)};
    if ((x[0] != x[1]) || (x[1] != x[2]) || (x[2] != x[3])) {
      fprintf(stderr, "Failed validation at %zd.\n", i);
      fprintf(stderr, "Simple:     %" PRIx64  "\n", x[0]);
      fprintf(stderr, "Recursive:  %" PRIx64  "\n", x[1]);
      fprintf(stderr, "Binary:     %" PRIx64  "\n", x[2]);
      fprintf(stderr, "Boosted<5>: %" PRIx64  "\n", x[3]);
      return 1;
    }
  }
  free(data);
  printf("Success up to %zd\n", length);
  return 0;
}
