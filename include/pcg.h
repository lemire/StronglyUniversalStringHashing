
#ifndef PCG_H_
#define PCG_H_



struct pcg_state_setseq_64 {    // Internals are *Private*.
    uint64_t state;             // RNG state.  All values are possible.
    uint64_t inc;               // Controls which RNG sequence (stream) is
    // selected. Must *always* be odd.
};
typedef struct pcg_state_setseq_64 pcg32_random_t;
typedef __uint128_t pcg128_t;
struct pcg_state_setseq_128 {
    pcg128_t state;
    pcg128_t inc;
};
typedef struct pcg_state_setseq_128 pcg64_random_t;

static pcg32_random_t pcg32_global = { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL };

#define PCG_128BIT_CONSTANT(high,low) \
        ((((pcg128_t)high) << 64) + low)

#define PCG_DEFAULT_MULTIPLIER_128 \
        PCG_128BIT_CONSTANT(2549297995355413924ULL,4865540595714422341ULL)

#define PCG64_INITIALIZER                                       \
    { PCG_128BIT_CONSTANT(0x979c9a98d8462005ULL, 0x7d3e9cb6cfe0549bULL),       \
      PCG_128BIT_CONSTANT(0x0000000000000001ULL, 0xda3e39cb94b95bdbULL) }

static pcg64_random_t pcg64_global = PCG64_INITIALIZER;

static inline uint32_t pcg32_random_r(pcg32_random_t* rng) {
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

static inline uint64_t pcg64_random_r(struct pcg_state_setseq_128* rng) {
    // the 32-bit version uses the old state to generate the next value
    rng->state = rng->state * PCG_DEFAULT_MULTIPLIER_128 + rng->inc;
    uint64_t xorshifted = ((uint64_t)(rng->state >> 64u)) ^ (uint64_t)rng->state;
    unsigned int rot = rng->state >> 122u;
    return (xorshifted >> rot) | (xorshifted << ((- rot) & 63));
}

static inline uint32_t pcg32_random() {
    return pcg32_random_r(&pcg32_global);
}

static inline uint64_t pcg64_random() {
    return pcg64_random_r(&pcg64_global);
}


#endif 
