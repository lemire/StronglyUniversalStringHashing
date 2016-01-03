#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include <unistd.h>

#include "simple-treehash.h"
#include "recursive-treehash.h"
#include "binary-treehash.h"


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
  ssize_t length = 1;
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
  uint64_t r[64];
  uint64_t * const data = malloc(sizeof(uint64_t) * length);
  if (!data) {
    fprintf(stderr, "Malloc failed: %zd\n", length);
    return 1;
  }
  if (!fill_random(r, 64) || !fill_random(data, length)) {
    return 1;
  }
  for (ssize_t i = 1; i < length; ++i) {
    const uint64_t x = simple_treehash(r, data, i),
      y = recursive_treehash(r, data, i),
      z = binary_treehash(r, data, i);
    if ((x != y) || (y != z)) {
      fprintf(stderr, "Failed validation at %zd.\n", i);
      fprintf(stderr, "Simple:    %0lx\n", x);
      fprintf(stderr, "Recursive: %0lx\n", y);
      fprintf(stderr, "Binary:    %0lx\n", z);
      return 1;
    }
  }
  puts("Success");
  return 0;
}
