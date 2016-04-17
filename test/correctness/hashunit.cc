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
#include "clmulhierarchical64bits.h"
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
  NAMED((&hashPMP64)),
  NAMED((&CLHASH)),
  NAMED((&hashCity)),
  NAMED((&hashSipHash)),
  NAMED((&hashVHASH64)),
  NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, CLNH, 7>)),
  NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, NHCL, 7>)),
  NAMED((&generic_treehash<BoostedZeroCopyGenericBinaryTreehash, NH, 7>))
};

const int HowManyFunctions64 =
    sizeof(hashFunctions) / sizeof(hashFunctions[0]);


// flipping a bit and then flipping it back should keep things the same
int testbitflipping() {
    printf("[%s] %s\n", __FILE__, __func__);
    int i = 0;
    int lengthStart = 1, lengthEnd = 1024; // inclusive
    uint64_t randbuffer[150] __attribute__ ((aligned (16)));// 150 should be plenty

    uint64_t * intstring;
    void * intstringoffsetted; // on purpose, we mess with the alignment
    intstringoffsetted = malloc(sizeof(uint64_t)*lengthEnd  + 1);
    if(intstringoffsetted == NULL) {
      cerr << "Failed to allocate " << lengthEnd << " words." << endl;
      return 1;
    }
    intstring =  (uint64_t *) ((char *) intstringoffsetted + 1);
    bool stillinplay [HowManyFunctions64];
    for (i = 0; i < HowManyFunctions64; ++i) {
        stillinplay[i] = true;
    }
    int re = 0;
    for(int trial = 0; trial < 50 ; ++ trial ) {
        for (i = 0; i < 150; ++i) {
            randbuffer[i] = pcg64_random();
        }
        for (i = 0; i < lengthEnd; ++i) {
            intstring[i] = pcg64_random();
        }
        for (i = 0; i < HowManyFunctions64; ++i) {
            if(! stillinplay[i]) continue;
            const hashFunction64 thisfunc64 = hashFunctions[i].f;
            for (int word = lengthStart; word < lengthEnd; ++word) {
                for (unsigned int bit = 0; bit  < sizeof(uint64_t); ++bit) {
                    const uint64_t first = thisfunc64(randbuffer, intstring, lengthEnd);
                    intstring[word] ^= (UINT64_C(1) << bit);
                    const uint64_t second = thisfunc64(randbuffer, intstring, lengthEnd);
                    intstring[word] ^= (UINT64_C(1) << bit);
                    const uint64_t third = thisfunc64(randbuffer, intstring, lengthEnd);
                    if (third != first) {
                        cout << "testing " << hashFunctions[i].name << endl;
                        cerr << "You have a bug." << endl;
                        cerr << endl;
                        stillinplay[i] = false;
                        re = 1;
                        goto endoflength;
                    }
                    if (second == first) {
                        cout << "testing " << hashFunctions[i].name << endl;
                        cerr << "You may have a bug. Flipping a bit did not change the value! " << endl;
                        cerr << endl;
                        stillinplay[i] = false;
                        re = 1;
                        goto endoflength;
                    }
                }
            }
endoflength:
            {}
        }
    }
    free(intstringoffsetted);
    return re;
}

int testunused() {
    printf("[%s] %s\n", __FILE__, __func__);
    assert(HowManyFunctions64 <= 64);
    int lengthStart = 1, lengthEnd = 1024; // inclusive
    int i;
    int length;
    uint64_t randbuffer[150] __attribute__ ((aligned (16)));// 150 should be plenty

    uint64_t * intstring = (uint64_t *) malloc(sizeof(uint64_t)*lengthEnd);
    for (i = 0; i < 150; ++i) {
        randbuffer[i] = pcg64_random();
    }
    for (i = 0; i < lengthEnd; ++i) {
        intstring[i] = pcg64_random();
    }
    int result = 0;
    for (i = 0; i < HowManyFunctions64; ++i) {
        const hashFunction64 thisfunc64 = hashFunctions[i].f;
        cout << "testing " << hashFunctions[i].name << endl;
        for (length = lengthStart; length <= lengthEnd; length += 1) {
            for (int place = 0; place < length; ++place) {
                const auto first_run = thisfunc64(randbuffer, intstring, length);
                const auto old_val = intstring[place];
                const uint64_t new_val = old_val + 1;
                intstring[place] = new_val;
                const auto second_run = thisfunc64(randbuffer, intstring, length);
                if (first_run == second_run) {
                    const auto old_val2 = intstring[place];
                    const uint64_t new_val2 = old_val2 + 1;
                    intstring[place] = new_val2;
                    const auto second_run2 = thisfunc64(randbuffer, intstring, length);
                    if (second_run2 == second_run) {
                        cerr << " You may have a bug. \n" <<endl;
                        cerr << second_run << " " << second_run2 << endl
                             << new_val << " " << new_val2 << endl;
                        cerr << endl;
                        result = 1;
                        goto endoflength;
                    }
                }
            }
        }
endoflength:
        {}
    }
    free(intstring);
    return result;
}

int main(int c, char ** arg) {
    (void) (c);
    (void) (arg);
    int r = 0;
    r |= testbitflipping();
    r |= testunused();
    if(r == 0) cout <<" Your code is probably ok." <<endl;
    else cout << "Your code looks buggy." << endl;
    return r;
}
