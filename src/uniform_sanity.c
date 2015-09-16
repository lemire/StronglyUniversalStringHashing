#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include "hashfunctions64bits.h"
#include "clmulhierarchical64bits.h"
#include "sys/time.h"

/** a quick sanity check for City, pclmul, vhash
 * to see that that a few randomly chosen members
 * of the family spread their outputs around
 * 1. A broken implementation still might do this
 * 2. A correct implementation might not (low probability)
 *
 * Assumption is that many bugs might result in uneven
 * distributions.
 **/

#define BITS_TO_BUCKET 8
#define NBUCKETS (1<< BITS_TO_BUCKET)
#define N 1000000
#define STRLEN 512
#define NTRIALS 50

void bucketize(uint64 hashval, int counts[], int positions[]) {
    int i, index = 0;

    for (i=0; i < BITS_TO_BUCKET; ++i) {
        int bit = ( hashval & (1L << positions[i])) ? 1 : 0;
        index = (index << 1) | bit;
    }
    counts[index]++;
}

uint64 rand64() {
    uint64 hi_bits =  ( (uint64) rand()) << 32;
    return rand() | hi_bits;
}

void check_distribution( int counts[]) {
    int min_count = N, max_count = 0;
    for (int i=0; i < NBUCKETS; ++i) {
        if (min_count > counts[i]) min_count = counts[i];
        if (max_count < counts[i]) {
            max_count = counts[i];
        }
    }

    if (min_count * 1.3 + 2 < max_count)
        printf("fishy count distribution %d %d\n", min_count, max_count);
}


int main() {
    int i,k;
    int bucket_bit_positions[BITS_TO_BUCKET];
    int bucket[NBUCKETS];
    uint64 *keys, *rands;
    hashFunction64 hfcns[3]= {&hashCity, &hashVHASH64, &CLHASH};
    struct timeval begin,end;

    keys  = calloc(N,sizeof(uint64));
    rands = calloc(N,sizeof(uint64)); // is it aligned for 128-bit fetch? spec: Supposed to be aligned for ANY data size, so ok...

    for (k=0; k < NTRIALS; ++k) {

        //uint64 small_step = 1 + (rand() % 100);
        //if (k==0) small_step = 1;

        for (i=0; i < N; ++i) {
            //printf ("%d\n",i);
            keys[i] = rand64(); // r+i*small_step;
            rands[i] = rand64();
        }

        // choose some (distinct)  random bit positions to inspect
        for (i=0; i < BITS_TO_BUCKET; ++i) {
            int bitpos, already_used;
            do {
                bitpos = rand() & 63;
                already_used=0;
                for (int ii=0; ii < i; ++ii)
                    if (bucket_bit_positions[ii] == bitpos) {
                        already_used = 1;
                        break;
                    }
            } while (already_used);
            bucket_bit_positions[i] = bitpos;
        }

        for (int algo=0; algo < 3; ++algo) {
            hashFunction64 fn = hfcns[algo];

            printf("** algo %d** ",algo);

            // init buckets
            for (i=0; i < NBUCKETS; ++i)
                bucket[i]=0;
            // lots of overlap between tested strings
            gettimeofday(&begin,NULL);

            for (int i=0; i < N-STRLEN; ++i)  {
                uint64 hval = fn(rands,  keys+i , STRLEN);  // STRLEN is in words, as it should be
                bucketize(hval, bucket, bucket_bit_positions);
            }
            gettimeofday(&end,NULL);
            printf("elapsed time is %g ms\n", ( (end.tv_sec - begin.tv_sec)*1e6+(end.tv_usec - begin.tv_usec))/1e3);
            check_distribution(bucket);
        }
    }
}
