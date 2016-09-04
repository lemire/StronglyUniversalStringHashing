#ifndef HASHFUNCTIONS32BITS_H_
#define HASHFUNCTIONS32BITS_H_

#include <immintrin.h>

// first pointer is a random source,
// next pointer is the data.
// outputs a hash value
typedef uint32_t (*hashFunction)(const void *  ,const  uint32_t * , const size_t );

//
// this is strongly universal. Condition: randomsource must be at least as long as
// the string length  + 3.
//
// Reference: Owen Kaser and Daniel Lemire, Strongly universal string hashing is fast, Computer Journal
// http://arxiv.org/abs/1202.4961
#if defined(__GNUC__) && !( defined(__clang__) || defined(__INTEL_COMPILER  ))
__attribute__((optimize("no-tree-vectorize")))
// GCC has buggy SSE2 code generation in some cases
// Thanks to Nathan Kurz for noticing that GCC 4.7 requires no-tree-vectorize to produce correct results.
#endif
uint32_t hashMultilinear(const void *  rs, const uint32_t *  string, const size_t length) {
    const uint32_t * const endstring = string + length;
    const uint64_t *  randomsource = (const uint64_t *) rs;
    uint64_t sum = *(randomsource++);
    for(; string!= endstring; ++randomsource,++string ) {
        sum+= (*randomsource *  (uint64_t)(*string)) ;
    }
    sum += *randomsource;
    return (int) (sum>>32);
}

#if defined(__GNUC__) && !( defined(__clang__) || defined(__INTEL_COMPILER  ))
__attribute__((optimize("no-tree-vectorize")))
// GCC has buggy SSE2 code generation in some cases
#endif
uint32_t hashMultilinear2by2(const void *  rs, const uint32_t *  string, const size_t length) {
    assert ( length / 2 * 2 == length );// length is even
    const uint64_t *  randomsource = (const uint64_t *) rs;
    const uint32_t * const endstring = string + length;
    uint64_t sum = *(randomsource++);
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        sum+= (*randomsource *  (uint64_t)(*string)) + (*(randomsource+1) *  (uint64_t)(*(string+1)))  ;
    }
    sum += *randomsource;
    return (int) (sum>>32);
}

#if defined(__GNUC__) && !( defined(__clang__) || defined(__INTEL_COMPILER  ))
__attribute__((optimize("no-tree-vectorize")))
// GCC has buggy SSE2 code generation in some cases
#endif
uint32_t hashMultilinearhalf(const void *  rs, const uint32_t *  string, const size_t length) {
    assert ( length / 2 * 2 == length );// length is even
    const uint32_t * const endstring = string + length;
    const uint64_t *  randomsource = (const uint64_t *) rs;
    uint64_t sum = *(randomsource++);
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        sum+= (*randomsource +  (uint64_t)(*string)) * (*(randomsource+1) +  (uint64_t)(*(string+1)))  ;
    }
    sum += *randomsource;
    return (int) (sum>>32);
}

#if defined(__GNUC__) && !( defined(__clang__) || defined(__INTEL_COMPILER  ))
__attribute__((optimize("no-tree-vectorize")))
// GCC has buggy SSE2 code generation in some cases
#endif
uint32_t hashMultilineardouble(const void *  rs, const uint32_t *  string, const size_t length) {
    assert ( length / 2 * 2 == length );// length is pair
    const uint32_t * const endstring = string + length;
    const uint64_t *  randomsource = (const uint64_t *) rs;
    uint64_t sum = *(randomsource++);
    uint64_t s2 = 0;
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        sum += *randomsource *  (uint64_t)(*string);
        s2 +=  *(randomsource+1) *  (uint64_t)(*(string+1))  ;
    }
    sum += *randomsource + s2;
    return (int) (sum>>32);
}



//Black, J.; Halevi, S.; Krawczyk, H.; Krovetz, T. (1999). "UMAC: Fast and Secure Message Authentication". Advances in Cryptology (CRYPTO '99)., Equation 1
// just high bits
uint32_t hashNH(const void *  randomsource, const uint32_t *  string, const size_t length) {
    assert ( length / 2 * 2 == length );// length is pair
    const uint32_t * const endstring = string + length;
    uint64_t sum = 0;
    const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    for(; string!= endstring; randomsource32+=2,string+=2 ) {
        sum+=
            (uint64_t) ( *randomsource32+ *string) *
            (*(randomsource32 + 1) + *(string+1));
    }
    return sum>>32;
}



// Fowler-Noll-Vo hashing
// L. C. Noll, Fowler/Noll/Vo Hash, http://www.isthe.com/chongo/tech/comp/fnv/
uint32_t hashFNV1a(const void * useless ,const uint32_t *  string, const size_t length) {
    (void) (useless);
    int sum = 0;
    const uint32_t someprime = 31;
    for(size_t i = 0; i < length; ++i) {
        sum = ( (string[i]) ^ sum) * someprime ;
    }
    return sum;
}

// M. V. Ramakrishna, J. Zobel, Performance in practice of string hashing functions,
//in: Proc. Int. Conf. on Database Systems for Advanced Applications, 1997.
uint32_t hashSAX(const void * useless,  const uint32_t *  string, const size_t length) {
    (void) (useless);
    int sum = 0;
    const int L = 3;
    const int R = 5;
    for(size_t i = 0; i < length; ++i) {
        sum = sum ^((sum<<L)+(sum>>R) + string[i]) ;
    }
    return sum;
}


// Rabin-Karp Hashing
uint32_t hashRabinKarp(const void * useless, const uint32_t *  string, const size_t length) {
    (void) (useless);
    int sum = 0;
    const int someprime = 31;
    for(size_t i = 0; i < length; ++i) {
        sum = someprime * sum +  string[i] ;
    }
    return sum;
}

// This is similar to Rabin-Karp, but we avoid the multiplication
//D. J. Bernstein, CDB-Constant Database, http://cr.yp.to/cdb.html
uint32_t hashBernstein(const void * useless,const uint32_t *  string, const size_t length) {
    (void) (useless);
    int sum = 0;
    const int L = 3;
    for(size_t i = 0; i < length; ++i) {
        sum= (( sum<< L) + sum) ^ string[i] ;
    }
    return sum;
}

// Fowler-Noll-Vo hashing
// L. C. Noll, Fowler/Noll/Vo Hash, http://www.isthe.com/chongo/tech/comp/fnv/
uint32_t hashFNV1(const void * useless,const uint32_t *  string, const size_t length) {
    (void) (useless);
    int sum = 0;
    const uint32_t someprime = 31;
    for(size_t i = 0; i < length; ++i) {
        sum = ( someprime * sum) ^ string[i] ;
    }
    return sum;
}














// used by pyramidal_Multilinear below
// code not thoroughly checked
void __hashMulti(const void *  randomsource, const uint32_t *  string, uint32_t *  output, int length, int blocksize) {
    int bufferpos = 0;
    int offset = length - length/blocksize*blocksize;
    for(; bufferpos<length/blocksize; ++bufferpos) {
        output[bufferpos] = hashMultilinear(randomsource, string+bufferpos*blocksize, blocksize);
    }
    if(offset>0) {
        output[length/blocksize] = hashMultilinear(randomsource, string+length/blocksize*blocksize,
                                   offset);
    }
}


// this function is 4/2**32 almost universal on 32 bits
// it uses at most 8KB of random bits (for strings having 32-bit lengths)
// code not thoroughly checked
uint32_t pyramidal_Multilinear(const void *  randomsource, const uint32_t *  string, const size_t len) {
    size_t length = len;
    int blocksize = 256;
    int newlength = ((length+blocksize-1)/blocksize);
    uint32_t * array = (uint32_t *)  malloc(newlength * sizeof(uint32_t));
    __hashMulti(randomsource, string, array,  length, blocksize);
    const char * randomChars = (const char *)randomsource;
    randomChars+=blocksize+1;
    length = newlength;
    while(length>1) {
        newlength = ((length+blocksize-1)/blocksize);
        uint32_t * array2 = (uint32_t *) malloc(newlength * sizeof(uint32_t));
        __hashMulti(randomChars, array, array2,  length, blocksize);
        randomChars+=blocksize+1;
        free(array);
        array = array2;
        length = newlength;
        if(length == 1) {
        }
    }
    uint32_t answer = array[0];
    free(array);
    return answer;
}

// Almost-strongly universal pseudo dot product (aka NH, aka half multilinear) using AVX
// and AVX2 instructions.
uint32_t pdp32avx(const void *rs, const uint32_t *string, const size_t length) {
  const uint32_t *const endstring = string + length;
  const __m256i *randomsource = (const __m256i *)rs;
  assert(((uintptr_t)randomsource & 31) == 0);
  assert ((length & 7) == 0);
  __m256i acc = *randomsource;
  ++randomsource;
  for (; string + 7 < endstring; randomsource += 1, string += 8) {
    __m256i input = _mm256_loadu_si256((const __m256i *)string);
    input = _mm256_add_epi32(input, *randomsource);
    __m256i hi = _mm256_srli_epi64(input, 32);
    input = _mm256_mul_epu32(input, hi);
    acc = _mm256_add_epi64(acc, input);
  }
  const int64_t intermediate[2] = {(int64_t)length, _mm256_extract_epi64(acc, 0) +
          _mm256_extract_epi64(acc, 1) + _mm256_extract_epi64(acc, 2) +
          _mm256_extract_epi64(acc, 3)};
  return hashMultilinear(randomsource, (const uint32_t *)intermediate, 4);
}

#endif /* HASHFUNCTIONS32BITS_H_ */
