
#ifndef MERSENNE_H_
#define MERSENNE_H_

#define NMERSENNE 624
#define MMERSENNE 397

typedef struct  {
    unsigned int MT[NMERSENNE + 1];
    unsigned int *map[NMERSENNE];
    int nValues;
} ZRandom;

inline static void initZRandom(ZRandom * ZR, unsigned iSeed)  {
    ZR->nValues = 0;
    // Seed the array used in random number generation.
    ZR->MT[0] = iSeed;
    for (int i = 1; i < NMERSENNE; ++i) {
        ZR->MT[i] = 1 + (69069 * ZR->MT[i - 1]);
    }
    // Compute map once to avoid % in inner loop.
    for (int i = 0; i < NMERSENNE; ++i) {
        ZR->map[i] = ZR->MT + ((i + MMERSENNE) % NMERSENNE);
    }
}



inline static unsigned int getValue(ZRandom * ZR) {
    if (0 == ZR->nValues) {
        ZR->MT[NMERSENNE] = ZR->MT[0];
        for (int i = 0; i < NMERSENNE; ++i) {
            unsigned y = (0x80000000 & ZR->MT[i]) | (0x7FFFFFFF
                         & ZR->MT[i + 1]);
            unsigned v = *(ZR->map[i]) ^ (y >> 1);
            if (1 & y)
                v ^= 2567483615;
            ZR->MT[i] = v;
        }
        ZR->nValues = NMERSENNE;
    }
    unsigned y = ZR->MT[NMERSENNE - ZR->nValues--];
    y ^= y >> 11;
    y ^= (unsigned int)((y << 7) & 2636928640);
    y ^= (unsigned int)((y << 15) & 4022730752);
    y ^= y >> 18;
    return y;
}


#endif /* MERSENNE_H_ */
