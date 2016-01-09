/////////////////////////////////////
// This C code is a companion to the paper
//
/////////////////////////////////////

//
// this code will hash strings of 64-bit characters. To use on
// strings of 8-bit characters, you may need some adequate padding.
//
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <assert.h>

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
#include "bigendianuniversal.h"
}

#include "treehash/simple-treehash.hh"
#include "treehash/recursive-treehash.hh"
#include "treehash/binary-treehash.hh"
#include "treehash/boosted-treehash.hh"

#define HowManyFunctions64 21

hashFunction64 funcArr64[HowManyFunctions64] = { &hashVHASH64, &CLHASH,
                                                 &hashCity, &hashSipHash,&GHASH64bit
                                                ,&hornerHash
                                                ,&unrolledHorner4
                                                ,&twiceHorner32
                                                 ,&iterateCL11
                                                 ,&treeCL9
                                                 ,&simple_treehash
                                                 ,&recursive_treehash
                                                 ,&binary_treehash
                                                 ,&boosted_treehash<1>
                                                 ,&boosted_treehash<2>
                                                 ,&boosted_treehash<3>
                                                 ,&boosted_treehash<4>
                                                 ,&boosted_treehash<5>
                                                 ,&boosted_treehash<6>
                                                 ,&boosted_treehash<7>
                                                 ,&simple_cl_treehash
                                               };

const char* functionnames64[HowManyFunctions64] = { "64-bit VHASH        ",
                                                    "64-bit CLHASH       ", "Google's City       ", "SipHash             ","GHASH          ",
                                                    "hornerHash          ",
                                                    "unrolled Horner     ",
                                                    "twice Horner32      "
                                                    ,"iterateCL 11        "
                                                    ,"treeCL9             "
                                                    ,"simple_treehash     "
                                                    ,"recursive_treehash  "
                                                    ,"binary_treehash     "
                                                    ,"boosted_treehash<1> "
                                                    ,"boosted_treehash<2> "
                                                    ,"boosted_treehash<3> "
                                                    ,"boosted_treehash<4> "
                                                    ,"boosted_treehash<5> "
                                                    ,"boosted_treehash<6> "
                                                    ,"boosted_treehash<7> "
                                                    ,"simple_cl_treehash  "

                                                  };

int main(int c, char ** arg) {
    (void) (c);
    (void) (arg);
    int which_algos = 0xffffffff;
    assert(HowManyFunctions64 <= 32);
    if (c > 1)
        which_algos = atoi(arg[1]);  // bitmask
    int lengthStart = 1, lengthEnd = 2048; // inclusive
    if (c > 2)
        lengthStart = atoi(arg[2]);
    if (c > 3)
        lengthEnd = atoi(arg[3]);

    int i, j;
    int length;
    int SHORTTRIALS;
    struct timeval start, finish;
    uint64_t randbuffer[150] __attribute__ ((aligned (16)));// 150 should be plenty
    uint32_t sumToFoolCompiler = 0;
    uint64_t * intstring;
    // We need 32 bytes of alignment for working with __m256i's
    if (posix_memalign((void **)(&intstring), 32, sizeof(uint64_t)*lengthEnd)) {
      cerr << "Failed to allocate " << lengthEnd << " words." << endl;
      return 1;
    }
    for (i = 0; i < 150; ++i) {
        randbuffer[i] = rand() | ((uint64_t)(rand()) << 32);
    }
    for (i = 0; i < lengthEnd; ++i) {
        intstring[i] = rand() | ((uint64_t)(rand()) << 32);
    }
    printf("#Reporting the number of cycles per byte.\n");
    printf("#First number is input length in  8-byte words.\n");
    printf("#          ");
    for (i = 0; i < HowManyFunctions64; ++i) {
        if (which_algos & (0x1 << i))
            printf("%s ", functionnames64[i]);
    }
    printf("\n");
    fflush(stdout);
    for (length = lengthStart; length <= lengthEnd; length += 1) {
        SHORTTRIALS = 8000000 / length;
        printf("%8d \t\t", length);
        hashFunction64 thisfunc64;
        for (i = 0; i < HowManyFunctions64; ++i) {
            if (!(which_algos & (0x1 << i)))
                continue;  // skip unselected algos
            thisfunc64 = funcArr64[i];
            sumToFoolCompiler += thisfunc64(randbuffer, intstring, length); // we do not count the first one
            gettimeofday(&start, 0);
            ticks lowest = ~(ticks)0;
            for (j = 0; j < SHORTTRIALS; ++j) {
                const ticks bef = startRDTSC();
                sumToFoolCompiler += thisfunc64(randbuffer, intstring, length);
                const ticks aft = stopRDTSCP();
                const ticks diff = aft-bef;
                lowest = (lowest < diff) ? lowest : diff;
            }
            gettimeofday(&finish, 0);
            printf(" %f ", lowest * 1.0 / (8.0 * length));
            fflush(stdout);
        }
        printf("\n");
    }
    free(intstring);
    printf("# ignore this #%d\n", sumToFoolCompiler);

}

