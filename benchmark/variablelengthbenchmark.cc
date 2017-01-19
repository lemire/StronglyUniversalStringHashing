/////////////////////////////////////
// This C code is a companion to the paper
//
/////////////////////////////////////

//
// this code will hash strings of 64-bit characters. To use on
// strings of 8-bit characters, you may need some adequate padding.
//
#include <cassert>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <sys/time.h>

#include <iostream>

using namespace std;

#ifdef __AVX__
#define __PCLMUL__ 1
#endif


extern "C" {
#include "timers.h"

#include "hashfunctions32bits.h"
#include "hashfunctions64bits.h"

#include "clmulhashfunctions32bits.h"
#include "clmulhashfunctions64bits.h"
#include "clmulpoly64bits.h"
#include "clmulhierarchical64bits.h"
#include "ghash.h"
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
    // From the 2015 paper:
    NAMED(&hashVHASH64), NAMED(&CLHASH), NAMED(&hashCity), NAMED(&hashSipHash),
    NAMED(&GHASH64bit),
    // Tree hashing:
    NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, CLNH, 7>)),
    NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, NHCL, 7>)),
    NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, NH, 7>)),
    NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, NHavx, 3>)),
    NAMED((&hashPMP64)),
};

const int HowManyFunctions64 =
    sizeof(hashFunctions) / sizeof(hashFunctions[0]);

int main(int c, char ** arg) {
    (void) (c);
    (void) (arg);
    bool aligned = true;
    if (NULL != strstr(arg[0], "unaligned")) {
      aligned = false;
    }
    uint64_t which_algos = ~0;
    assert(HowManyFunctions64 <= 64);
    if (c > 1) {
      if (1 != sscanf(arg[1], "%" SCNu64, &which_algos)) {
        return 1;
      }
    }
    int lengthStart = 1, lengthEnd = 2048; // inclusive
    if (c > 2)
        lengthStart = atoi(arg[2]);
    if (c > 3)
        lengthEnd = atoi(arg[3]);

    int i, j;
    int length;
    int SHORTTRIALS;
    struct timeval start, finish;
    uint64_t randbuffer[150] __attribute__ ((aligned (32)));// 150 should be plenty
    uint32_t sumToFoolCompiler = 0;
    uint64_t * intstring;
    // We need 32 bytes of alignment for working with __m256i's
    ++lengthEnd;
    if (posix_memalign((void **)(&intstring), 32, sizeof(uint64_t)*lengthEnd)) {
      cerr << "Failed to allocate " << lengthEnd << " words." << endl;
      return 1;
    }
    --lengthEnd;
    if (!aligned) {
      intstring = reinterpret_cast<decltype(intstring)>(1 + reinterpret_cast<char *>(intstring));
    }
    for (i = 0; i < 150; ++i) {
        randbuffer[i] = rand() | ((uint64_t)(rand()) << 32);
    }
    for (i = 0; i < lengthEnd; ++i) {
        intstring[i] = rand() | ((uint64_t)(rand()) << 32);
    }
    printf("#Reporting the number of cycles per byte.\n");
    printf("#First number is input length in  8-byte words.\n");
    printf("0 ");
    for (i = 0; i < HowManyFunctions64; ++i) {
        if (which_algos & (0x1ull << i))
          cout << '"' << hashFunctions[i].name << "\" ";
    }
    printf("\n");
    fflush(stdout);
    for (length = lengthStart; length <= lengthEnd; length += 1) {
        SHORTTRIALS = 8000000 / length;
        printf("%8d \t\t", length);

        for (i = 0; i < HowManyFunctions64; ++i) {
            if (!(which_algos & (0x1ull << i)))
                continue;  // skip unselected algos
            const hashFunction64 thisfunc64 = hashFunctions[i].f;
            sumToFoolCompiler += thisfunc64(randbuffer, intstring, length); // we do not count the first one
            gettimeofday(&start, 0);
            const ticks bef = startRDTSC();
            for (j = 0; j < SHORTTRIALS; ++j) {
                sumToFoolCompiler += thisfunc64(randbuffer, intstring, length);
            }
            const ticks aft = stopRDTSCP();
            gettimeofday(&finish, 0);
            printf(" %.3f ", ((aft-bef) * 1.0)/(8.0 * SHORTTRIALS * length));
            fflush(stdout);
        }
        printf("\n");
    }
    if (!aligned) {
      intstring = reinterpret_cast<decltype(intstring)>(-1 + reinterpret_cast<char *>(intstring));
    }
    free(intstring);
    printf("# ignore this #%d\n", sumToFoolCompiler);

}
