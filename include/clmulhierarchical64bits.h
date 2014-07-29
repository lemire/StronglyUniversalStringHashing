/*
 * clmulhierarchical64bits.h
 *
 *  Created on: Jul 29, 2014
 *      Author: lemire
 */

#ifndef CLMULHIERARCHICAL64BITS_H_
#define CLMULHIERARCHICAL64BITS_H_


// simplified version of hashGaloisFieldfast64halfunrolled_precomp for use with hashCLMULHierarchical
uint64_t __clmulhalfscalarproduct(const void*  rs, const uint64_t *  string, const size_t length) {
    assert(length / 4 * 4 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )rs;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string+3 < endstring; randomsource+=4,string+=4 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        const __m128i temp12 = _mm_lddqu_si128((__m128i * )(randomsource + 2));
        const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string+2));
        const __m128i add12 =  _mm_xor_si128 (temp12,temp22);
        const __m128i clprod12  = _mm_clmulepi64_si128( add12, add12, 0x10);
        acc = _mm_xor_si128 (clprod12,acc);
    }
    return precompReduction64(acc);
}

// simplified version of hashGaloisFieldfast64halfunrolled_precomp for use with hashCLMULHierarchical
uint64_t __clmulhalfscalarproductwithtail(const void*  rs, const uint64_t *  string, const size_t length) {
    const uint64_t * const endstring = string + length;
    const uint64_t *  randomsource = ( const uint64_t * )rs;
    __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string+3 < endstring; randomsource+=4,string+=4 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        const __m128i temp12 = _mm_lddqu_si128((__m128i * )(randomsource + 2));
        const __m128i temp22 = _mm_lddqu_si128((__m128i *) (string+2));
        const __m128i add12 =  _mm_xor_si128 (temp12,temp22);
        const __m128i clprod12  = _mm_clmulepi64_si128( add12, add12, 0x10);
        acc = _mm_xor_si128 (clprod12,acc);
    }
    if(string+1< endstring) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);
        randomsource+=2;
        string+=2;
    }
    if(string < endstring) {
        const __m128i temp1 = _mm_set_epi64x(0,*randomsource);
        const __m128i temp2 = _mm_set_epi64x(0,*string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        acc = _mm_xor_si128 (clprod1,acc);
    }
    return precompReduction64(acc);
}

// This should hash arbitrary strings fairly fast using no more than (256+1)*8 keys.
uint64_t hashCLMULHierarchical256(const void* rs, const uint64_t * string,
		const size_t length) {
        if(length == 0) return 0; // hmmmm...
	const int m = 256;
	const int l = 7;
	// trick here is that m^(l+1) is 2^64 which ensures we can support strings up to 2^64/8 easily.
	uint64_t hashtree[l][m];
	int counters[l];
	memset(counters,0,l*sizeof(int));
	size_t t = 0;
	for(; t + m <=  length; t+= m){
		hashtree[0][++counters[0]] = __clmulhalfscalarproduct(rs,string+t,m);
		// next, we push this up the tree, if needed
		for(int j = 0; counters[j] == m; ++j) {
			counters[j] = 0;
			hashtree[j+1][++counters[j+1]] = __clmulhalfscalarproduct(rs+(m+1)*(j+1),hashtree[j],m);
		}
	}
	int leftover = length - t;
	if(leftover > 0) {
		hashtree[0][++counters[0]] = __clmulhalfscalarproductwithtail(rs,string+t,leftover);
	}
	// next, we push this up the tree, if needed
        int maxj = 0;
        for(int j = 0; j < l ; ++j) if(counters[j] > 0) maxj = j;
	for(int j = 0; j < maxj; ++j) {
		if(counters[j] > 0) {
			 hashtree[j+1][++counters[j+1]] = __clmulhalfscalarproductwithtail(rs+(m+1)*(j+1),hashtree[j],counters[j]);
			 counters[j] = 0;
		}
	}
        if(counters[maxj] == 1) return hashtree[maxj][0];
        return __clmulhalfscalarproductwithtail(rs+(m+1)*(maxj+1),hashtree[maxj],counters[maxj]);
}

// uses no more than (128+1) * 9 keys
uint64_t hashCLMULHierarchical128(const void* rs, const uint64_t * string,
		const size_t length) {
        if(length == 0) return 0; // hmmmm...
	const int m = 128;
	const int l = 8;
	// trick here is that m^(l+1) is 2^64 which ensures we can support strings up to 2^64/8 easily.
	uint64_t hashtree[l][m];
	int counters[l];
	memset(counters,0,l*sizeof(int));
	size_t t = 0;
	for(; t + m <=  length; t+= m){
		hashtree[0][++counters[0]] = __clmulhalfscalarproduct(rs,string+t,m);
		// next, we push this up the tree, if needed
		for(int j = 0; counters[j] == m; ++j) {
			counters[j] = 0;
			hashtree[j+1][++counters[j+1]] = __clmulhalfscalarproduct(rs+(m+1)*(j+1),hashtree[j],m);
		}
	}
	int leftover = length - t;
	if(leftover > 0) {
		hashtree[0][++counters[0]] = __clmulhalfscalarproductwithtail(rs,string+t,leftover);
	}
	// next, we push this up the tree, if needed
        int maxj = 0;
        for(int j = 0; j < l ; ++j) if(counters[j] > 0) maxj = j;
	for(int j = 0; j < maxj; ++j) {
		if(counters[j] > 0) {
			 hashtree[j+1][++counters[j+1]] = __clmulhalfscalarproductwithtail(rs+(m+1)*(j+1),hashtree[j],counters[j]);
			 counters[j] = 0;
		}
	}
        if(counters[maxj] == 1) return hashtree[maxj][0];
        return __clmulhalfscalarproductwithtail(rs+(m+1)*(maxj+1),hashtree[maxj],counters[maxj]);
}




#endif /* CLMULHIERARCHICAL64BITS_H_ */
