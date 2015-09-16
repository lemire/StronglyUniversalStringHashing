/*
 * These are polynomial hashing functions.
 */

#ifndef CLMULPOLY64BITS_H_
#define CLMULPOLY64BITS_H_

#include "clmul.h"


// given one 8-byte key, generate 128 bytes of powers (randomp, randomp**2, ...) or 2 cache lines.
// could be used in conjunction with CLMULPoly64CL2
void precomputePowers(uint64_t randomp, __m128i * out) {
    assert(((uintptr_t) out & 15) == 0); // we expect cache line alignment for the keys
    __m128i power1 = _mm_loadl_epi64((__m128i  const *) &randomp );
    __m128i power2 = precompReduction64_si128(
                         _mm_clmulepi64_si128(power1, power1, 0x00));
    // save power1, power2 in the first 128 bits
    _mm_store_si128(out++,_mm_unpacklo_epi64(power1,power2));
    __m128i powerodd = power1;
    __m128i powereven = power2;
    for(int k = 0; k<7; ++k) {
        // we use the fact that if we have p**k1, p**(k1+1), we can get p**(k1+2), p**(k1+3)
        // by multiply both by p**2
        powerodd =  precompReduction64_si128(
                        _mm_clmulepi64_si128(powerodd, power2, 0x00));
        powereven =  precompReduction64_si128(
                         _mm_clmulepi64_si128(powereven, power2, 0x00));
        _mm_store_si128(out++,_mm_unpacklo_epi64(powerodd,powereven));
    }
}





//CLMULPoly64CL2
uint64_t CLMULPoly64CL2(const void* rs, const uint64_t * string,
                        const size_t length) {
    const __m128i * randomsource = (const __m128i *) rs;
    // 4 128-bit is a cache line!!!
    assert(((uintptr_t) rs & 15) == 0); // we expect cache line alignment for the keys
    if (length == 0)
        return 0; // duh!

//	if p is some 64-bit key, the numbered keys below should contain
//	its reduced powers as follows
//	key1=(p,p**2)
//	key2=(p**3,p**4)
//	...
//	and so on

    const __m128i key1 = _mm_load_si128(randomsource);
    const __m128i key2 = _mm_load_si128(randomsource + 1);
    const __m128i key3 = _mm_load_si128(randomsource + 2);
    const __m128i key4 = _mm_load_si128(randomsource + 3);

    const uint64_t * const endstring = string + length;
    __m128i acc;
    if((length & 1) == 0) {// even
        acc= _mm_setzero_si128();
    } else {
        acc= _mm_loadl_epi64((__m128i  const *) string);
        ++string;
    }
    if (string + 15 < endstring) {
        const __m128i key5 = _mm_load_si128(randomsource + 4);
        const __m128i key6 = _mm_load_si128(randomsource + 5);
        const __m128i key7 = _mm_load_si128(randomsource + 6);
        const __m128i key8 = _mm_load_si128(randomsource + 7);

        for (; string + 15 < endstring; string += 16) {
            __m128i p1 = _mm_clmulepi64_si128(acc, key8, 0x00);

            __m128i temp1 = _mm_lddqu_si128((__m128i *) string);
            __m128i q1 = _mm_xor_si128(temp1, key1);
            __m128i p2 = _mm_clmulepi64_si128(q1, q1, 0x01);

            __m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2));
            __m128i q2 = _mm_xor_si128(temp2, key2);
            __m128i p3 = _mm_clmulepi64_si128(q2, q2, 0x01);

            __m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4));
            __m128i q3 = _mm_xor_si128(temp3, key3);
            __m128i p4 = _mm_clmulepi64_si128(q3, q3, 0x01);

            __m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6));
            __m128i q4 = _mm_xor_si128(temp4, key4);
            __m128i p5 = _mm_clmulepi64_si128(q4, q4, 0x01);

            __m128i temp5 = _mm_lddqu_si128((__m128i *) (string + 8));
            __m128i q5 = _mm_xor_si128(temp5, key5);
            __m128i p6 = _mm_clmulepi64_si128(q5, q5, 0x01);

            __m128i temp6 = _mm_lddqu_si128((__m128i *) (string + 10));
            __m128i q6 = _mm_xor_si128(temp6, key6);
            __m128i p7 = _mm_clmulepi64_si128(q6, q6, 0x01);

            __m128i temp7 = _mm_lddqu_si128((__m128i *) (string + 12));
            __m128i q7 = _mm_xor_si128(temp7, key7);
            __m128i p8 = _mm_clmulepi64_si128(q7, q7, 0x01);

            __m128i temp8 = _mm_lddqu_si128((__m128i *) (string + 14));
            __m128i p9 = _mm_clmulepi64_si128(temp8, key8, 0x01);
            __m128i p10 = _mm_slli_si128(temp8, 8);

            __m128i Q1 = _mm_xor_si128(p1, p2);
            __m128i Q2 = _mm_xor_si128(p3, p4);
            __m128i Q3 = _mm_xor_si128(p5, p6);
            __m128i Q4 = _mm_xor_si128(p7, p8);
            __m128i Q5 = _mm_xor_si128(p9, p10);

            __m128i r1 = _mm_xor_si128(Q1, Q2);
            __m128i r2 = _mm_xor_si128(Q3, Q4);
            __m128i R1 = _mm_xor_si128(r1, r2);
            R1 = _mm_xor_si128(R1, Q5);
            acc = precompReduction64_si128(R1);
        }
    }
    if (string + 7 < endstring) {
        __m128i p1 = _mm_clmulepi64_si128(acc, key4, 0x00);

        __m128i temp1 = _mm_lddqu_si128((__m128i *) string);
        __m128i q1 = _mm_xor_si128(temp1, key1);
        __m128i p2 = _mm_clmulepi64_si128(q1, q1, 0x01);

        __m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2));
        __m128i q2 = _mm_xor_si128(temp2, key2);
        __m128i p3 = _mm_clmulepi64_si128(q2, q2, 0x01);

        __m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4));
        __m128i q3 = _mm_xor_si128(temp3, key3);
        __m128i p4 = _mm_clmulepi64_si128(q3, q3, 0x01);

        __m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6));
        __m128i p5 = _mm_clmulepi64_si128(temp4, key4, 0x01);
        __m128i p6 = _mm_slli_si128(temp4, 8);

        __m128i Q1 = _mm_xor_si128(p1, p2);
        __m128i Q2 = _mm_xor_si128(p3, p4);
        __m128i Q3 = _mm_xor_si128(p5, p6);

        __m128i r1 = _mm_xor_si128(Q1, Q2);
        __m128i r2 = _mm_xor_si128(Q3, r1);
        acc = precompReduction64_si128(r2);
        string += 8;
    }
    if (string + 3 < endstring) {
        __m128i p1 = _mm_clmulepi64_si128(acc, key2, 0x00);

        __m128i temp1 = _mm_lddqu_si128((__m128i *) string);
        __m128i q1 = _mm_xor_si128(temp1, key1);
        __m128i p2 = _mm_clmulepi64_si128(q1, q1, 0x01);

        __m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2));
        __m128i p3 = _mm_clmulepi64_si128(temp2, key2, 0x01);
        __m128i p4 = _mm_slli_si128(temp2, 8);

        __m128i Q1 = _mm_xor_si128(p1, p2);
        __m128i Q2 = _mm_xor_si128(p3, p4);

        __m128i r1 = _mm_xor_si128(Q1, Q2);
        acc = precompReduction64_si128(r1);
        string += 4;
    }
    if(string + 1 < endstring) {
        __m128i p1 = _mm_clmulepi64_si128(acc, key1, 0x00);
        __m128i temp1 = _mm_lddqu_si128((__m128i *) string);
        __m128i p2 = _mm_clmulepi64_si128(temp1, key1, 0x01);
        __m128i p3 = _mm_slli_si128(temp1, 8);
        acc = precompReduction64_si128(
                  _mm_xor_si128(_mm_xor_si128(p1, p2), p3));
        string += 2;
    }
    return _mm_cvtsi128_si64(acc);
}


#endif /* CLMULPOLY64BITS_H_ */

// Rest is crap to be deleted eventually

//CLMULPoly64CL1
//uint64_t CLMULPoly64CL1(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const __m128i * randomsource = (const __m128i *) rs;
//	// 4 128-bit is a cache line!!!
//	__m128i key1 = _mm_lddqu_si128( randomsource);
//	__m128i key2 = _mm_lddqu_si128( randomsource + 1 );
//	__m128i key3 = _mm_lddqu_si128( randomsource + 2 );
//	__m128i key4 = _mm_lddqu_si128( randomsource + 3 );
//
//	const uint64_t * const endstring = string + length;
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	for(;string + 7 < endstring;string+=8) {
//		__m128i p1 = _mm_clmulepi64_si128(acc, key4, 0x00);
//
//		__m128i temp1 = _mm_lddqu_si128((__m128i *) string);
//		__m128i q1 = _mm_xor_si128(temp1,key1);
//		__m128i p2 = _mm_clmulepi64_si128(q1, q1, 0x01);
//
//
//		__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2));
//		__m128i q2 = _mm_xor_si128(temp2,key2);
//		__m128i p3 = _mm_clmulepi64_si128(q2, q2, 0x01);
//
//		__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4));
//		__m128i q3 = _mm_xor_si128(temp3,key3);
//		__m128i p4 = _mm_clmulepi64_si128(q3, q3, 0x01);
//
//
//		__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6));
//		__m128i p5 = _mm_clmulepi64_si128(temp4, key4, 0x01);
//		__m128i p6 = _mm_slli_si128(temp4, 8);
//
//
//		__m128i Q1 = _mm_xor_si128(p1,p2);
//		__m128i Q2 = _mm_xor_si128(p3,p4);
//		__m128i Q3 = _mm_xor_si128(p5,p6);
//
//		__m128i r1 = _mm_xor_si128(Q1,Q2);
//		__m128i r2 = _mm_xor_si128(Q3,r1);
//		acc = precompReduction64_si128(r2);
//	}
//
//	for (;string + 1 < endstring;string+=2) {
//		__m128i p1 = _mm_clmulepi64_si128(acc, key4, 0x00);
//		__m128i temp1 = _mm_lddqu_si128((__m128i *) string);
//		__m128i p2 = _mm_clmulepi64_si128(temp1, key4, 0x01);
//		__m128i p3 = _mm_slli_si128(temp1, 8);
//		acc = precompReduction64_si128(_mm_xor_si128(_mm_xor_si128(p1,p2),p3));
//	}
//	if(string<endstring){
//		__m128i p1 = _mm_clmulepi64_si128(acc, key1, 0x00);
//		acc = precompReduction64_si128(_mm_xor_si128(p1,_mm_set_epi64x(0, *string)));
//	}
//	return _mm_cvtsi128_si64(acc);
//}

//
//// simple 64-bit polynomial hashing, uses only one key
//// not expected to be fast!
//uint64_t hashGaloisFieldPoly64(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const uint64_t * randomsource = (const uint64_t *) rs;
//	assert(*randomsource != 0); //otherwise silly
//	const uint64_t * const endstring = string + length;
//	__m128i key = _mm_set_epi64x(0, *(randomsource));
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	for (; string < endstring; ++string) {
//		const __m128i temp = _mm_set_epi64x(0, *string);
//		const __m128i multi = _mm_clmulepi64_si128(acc, key, 0x00);
//		acc = barrettWithoutPrecomputation64_si128(multi);
//		acc = _mm_xor_si128(acc, temp);
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
//uint64_t precomphashGaloisFieldPoly64(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const uint64_t * randomsource = (const uint64_t *) rs;
//	assert(*randomsource != 0); //otherwise silly
//	const uint64_t * const endstring = string + length;
//	__m128i key = _mm_set_epi64x(0, *(randomsource));
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	for (; string < endstring; ++string) {
//		const __m128i temp = _mm_set_epi64x(0, *string);
//		const __m128i multi = _mm_clmulepi64_si128(acc, key, 0x00);
//		acc = precompReduction64_si128(multi);
//		acc = _mm_xor_si128(acc, temp);
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
//uint64_t fasthashGaloisFieldPoly64_2_noprecomp(const void* rs,
//		const uint64_t * string, const size_t length) {
//	const uint64_t * randomsource = (const uint64_t *) rs;
//	assert(*randomsource != 0); //otherwise silly
//	const uint64_t * const endstring = string + length;
//	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	if (string + 1 < endstring) {
//		// we start by precomputing the powers of the key
//		__m128i tkey2 = barrettWithoutPrecomputation64_si128(
//				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
//		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
//				_mm_slli_si128(tkey2, 8));
//		for (; string + 1 < endstring; string += 2) {
//			__m128i temp = _mm_lddqu_si128((__m128i *) string);
//			acc = _mm_clmulepi64_si128(acc, key, 0x10);
//			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
//			acc = _mm_xor_si128(clprod1, acc);
//			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
//			acc = barrettWithoutPrecomputation64_si128(acc);
//		}
//	}
//	if (string < endstring) {
//		const __m128i temp = _mm_set_epi64x(0, *string);
//		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
//		acc = barrettWithoutPrecomputation64_si128(multi);
//		acc = _mm_xor_si128(acc, temp);
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
//uint64_t fasthashGaloisFieldPoly64_2(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const uint64_t * randomsource = (const uint64_t *) rs;
//	assert(*randomsource != 0); //otherwise silly
//	const uint64_t * const endstring = string + length;
//	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	if (string + 1 < endstring) {
//		// we start by precomputing the powers of the key
//		__m128i tkey2 = precompReduction64_si128(
//				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
//		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
//				_mm_slli_si128(tkey2, 8));
//		for (; string + 1 < endstring; string += 2) {
//			__m128i temp = _mm_lddqu_si128((__m128i *) string);
//			acc = _mm_clmulepi64_si128(acc, key, 0x10);
//			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
//			acc = _mm_xor_si128(clprod1, acc);
//			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
//			acc = precompReduction64_si128(acc);
//		}
//	}
//	if (string < endstring) {
//		const __m128i temp = _mm_set_epi64x(0, *string);
//		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
//		acc = precompReduction64_si128(multi);
//		acc = _mm_xor_si128(acc, temp);
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
//// fast 64-bit polynomial hashing, uses only one key
//// expected to be fast!
////TODO: can use more keys for increased universality
//uint64_t fasthashGaloisFieldPoly64_4(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const uint64_t * randomsource = (const uint64_t *) rs;
//	assert(*randomsource != 0); //otherwise silly
//	const uint64_t * const endstring = string + length;
//	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	if (string + 1 < endstring) {
//		// we start by precomputing the powers of the key
//		__m128i tkey2 = precompReduction64_si128(
//				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
//		// powers of the keys are packed into two registers
//		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
//				_mm_slli_si128(tkey2, 8));
//		if (string + 3 < endstring) {
//			__m128i tkey3 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
//			__m128i tkey4 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
//			__m128i key2 = _mm_xor_si128(
//					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
//					_mm_slli_si128(tkey4, 8));
//			for (; string + 3 < endstring; string += 4) {
//				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
//				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
//				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
//				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
//				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
//				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
//				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
//				acc = precompReduction64_si128(
//						_mm_xor_si128(acc,
//								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
//										_mm_xor_si128(clprod2, clprod3))));
//			}
//		}
//		for (; string + 1 < endstring; string += 2) {
//			__m128i temp = _mm_lddqu_si128((__m128i *) string);
//			acc = _mm_clmulepi64_si128(acc, key, 0x10);
//			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
//			acc = _mm_xor_si128(clprod1, acc);
//			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
//			acc = precompReduction64_si128(acc);
//		}
//	}
//	if (string < endstring) {
//		const __m128i temp = _mm_set_epi64x(0, *string);
//		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
//		acc = precompReduction64_si128(multi);
//		acc = _mm_xor_si128(acc, temp);
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
//
//// fast 64-bit polynomial hashing, uses only one key
//// expected to be fast!
//uint64_t fasthashGaloisFieldPoly64_8(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const uint64_t * randomsource = (const uint64_t *) rs;
//	assert(*randomsource != 0); //otherwise silly
//	const uint64_t * const endstring = string + length;
//	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	if (string + 1 < endstring) {
//		// we start by precomputing the powers of the key
//		__m128i tkey2 = precompReduction64_si128(
//				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
//		// powers of the keys are packed into two registers
//		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
//				_mm_slli_si128(tkey2, 8));
//		if (string + 3 < endstring) {
//			__m128i tkey3 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
//			__m128i tkey4 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
//			__m128i key2 = _mm_xor_si128(
//					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
//					_mm_slli_si128(tkey4, 8));
//			if (string + 7 < endstring) {
//				__m128i tkey5 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey2, tkey3, 0x00));
//				__m128i tkey6 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey3, tkey3, 0x00));
//				__m128i tkey7 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey3, tkey4, 0x00));
//				__m128i tkey8 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey4, tkey4, 0x00));
//				__m128i key3 = _mm_xor_si128(
//						_mm_and_si128(tkey5, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey6, 8));
//				__m128i key4 = _mm_xor_si128(
//						_mm_and_si128(tkey7, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey8, 8));
//				for (; string + 7 < endstring; string += 8) {
//					__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
//					__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
//					__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4)); //a5 a6
//					__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6)); //a7 a8
//					const __m128i x1 = _mm_srli_si128(temp4, 8); //a8
//					const __m128i clprod1 = _mm_clmulepi64_si128(temp4, key,
//							0x00); //a7*k
//					const __m128i clprod2 = _mm_clmulepi64_si128(temp3, key,
//							0x11); //a6*k^2
//					const __m128i clprod3 = _mm_clmulepi64_si128(temp3, key2,
//							0x00); //a5*k^3
//					const __m128i clprod4 = _mm_clmulepi64_si128(temp2, key2,
//							0x11); //a4*k^4
//					const __m128i clprod5 = _mm_clmulepi64_si128(temp2, key3,
//							0x00); //a3*k^5
//					const __m128i clprod6 = _mm_clmulepi64_si128(temp, key3,
//							0x11); //a2*k^6
//					const __m128i clprod7 = _mm_clmulepi64_si128(temp, key4,
//							0x00); //a1*k^7
//					acc = _mm_clmulepi64_si128(acc, key4, 0x10); //k^8
//
//					const __m128i t1 = _mm_xor_si128(x1, clprod1);
//					const __m128i t2 = _mm_xor_si128(clprod2, clprod3);
//					const __m128i t3 = _mm_xor_si128(clprod4, clprod5);
//					const __m128i t4 = _mm_xor_si128(clprod6, clprod7);
//
//					const __m128i z1 = _mm_xor_si128(t1, t2);
//					const __m128i z2 = _mm_xor_si128(t3, t4);
//
//					acc = precompReduction64_si128(
//							_mm_xor_si128(acc, _mm_xor_si128(z1, z2)));
//				}
//
//			}
//			for (; string + 3 < endstring; string += 4) {
//				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
//				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
//				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
//				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
//				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
//				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
//				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
//				acc = precompReduction64_si128(
//						_mm_xor_si128(acc,
//								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
//										_mm_xor_si128(clprod2, clprod3))));
//			}
//		}
//		for (; string + 1 < endstring; string += 2) {
//			__m128i temp = _mm_lddqu_si128((__m128i *) string);
//			acc = _mm_clmulepi64_si128(acc, key, 0x10);
//			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
//			acc = _mm_xor_si128(clprod1, acc);
//			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
//			acc = precompReduction64_si128(acc);
//		}
//	}
//	if (string < endstring) {
//		const __m128i temp = _mm_set_epi64x(0, *string);
//		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
//		acc = precompReduction64_si128(multi);
//		acc = _mm_xor_si128(acc, temp);
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
//
//// experimental
//uint64_t fasthashGaloisFieldPoly64_16(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const uint64_t * randomsource = (const uint64_t *) rs;
//	assert(*randomsource != 0); //otherwise silly
//	const uint64_t * const endstring = string + length;
//	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	if (string + 1 < endstring) {
//		// we start by precomputing the powers of the key
//		__m128i tkey2 = precompReduction64_si128(
//				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
//		// powers of the keys are packed into two registers
//		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
//				_mm_slli_si128(tkey2, 8));
//		if (string + 3 < endstring) {
//			__m128i tkey3 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
//			__m128i tkey4 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
//			__m128i key2 = _mm_xor_si128(
//					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
//					_mm_slli_si128(tkey4, 8));
//			if (string + 15 < endstring) {
//				__m128i tkey5 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey2, tkey3, 0x00));
//				__m128i tkey6 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey3, tkey3, 0x00));
//				__m128i tkey7 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey3, tkey4, 0x00));
//				__m128i tkey8 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey4, tkey4, 0x00));
//				__m128i tkey9 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey5, tkey4, 0x00));
//				__m128i tkey10 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey5, tkey5, 0x00));
//				__m128i tkey11 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey6, tkey5, 0x00));
//				__m128i tkey12 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey6, tkey6, 0x00));
//				__m128i tkey13 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey6, tkey7, 0x00));
//				__m128i tkey14 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey7, tkey7, 0x00));
//				__m128i tkey15 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey8, tkey7, 0x00));
//				__m128i tkey16 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey8, tkey8, 0x00));
//				__m128i key3 = _mm_xor_si128(
//						_mm_and_si128(tkey5, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey6, 8));
//				__m128i key4 = _mm_xor_si128(
//						_mm_and_si128(tkey7, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey8, 8));
//				__m128i key5 = _mm_xor_si128(
//						_mm_and_si128(tkey9, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey10, 8));
//				__m128i key6 = _mm_xor_si128(
//						_mm_and_si128(tkey11, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey12, 8));
//				__m128i key7 = _mm_xor_si128(
//						_mm_and_si128(tkey13, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey14, 8));
//				__m128i key8 = _mm_xor_si128(
//						_mm_and_si128(tkey15, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey16, 8));
//
//				for (; string + 15 < endstring; string += 16) {
//					__m128i temp = _mm_lddqu_si128((__m128i *) string);
//					__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2));
//					__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4));
//					__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6));
//					__m128i temp5 = _mm_lddqu_si128((__m128i *) (string + 8));
//					__m128i temp6 = _mm_lddqu_si128((__m128i *) (string + 10));
//					__m128i temp7 = _mm_lddqu_si128((__m128i *) (string + 12));
//					__m128i temp8 = _mm_lddqu_si128((__m128i *) (string + 14));
//					 __m128i x1 = _mm_srli_si128(temp8, 8);
//					 __m128i clprod1 = _mm_clmulepi64_si128(temp8, key,
//							0x00);
//					 __m128i clprod2 = _mm_clmulepi64_si128(temp7, key,
//							0x11);
//					 __m128i clprod3 = _mm_clmulepi64_si128(temp7, key2,
//							0x00);
//					 __m128i clprod4 = _mm_clmulepi64_si128(temp6, key2,
//							0x11);
//					 __m128i clprod5 = _mm_clmulepi64_si128(temp6, key3,
//							0x00);
//					 __m128i clprod6 = _mm_clmulepi64_si128(temp5, key3,
//							0x11);
//					 __m128i clprod7 = _mm_clmulepi64_si128(temp5, key4,
//							0x00);
//					 __m128i clprod8 = _mm_clmulepi64_si128(temp4, key4,
//							0x11);
//					 __m128i clprod9 = _mm_clmulepi64_si128(temp4, key5,
//							0x00);
//					 __m128i clprod10 = _mm_clmulepi64_si128(temp3, key5,
//							0x11);
//					 __m128i clprod11 = _mm_clmulepi64_si128(temp3, key6,
//							0x00);
//					 __m128i clprod12 = _mm_clmulepi64_si128(temp2, key6,
//							0x11);
//					 __m128i clprod13 = _mm_clmulepi64_si128(temp2, key7,
//							0x00);
//					 __m128i clprod14 = _mm_clmulepi64_si128(temp, key7,
//							0x11);
//					 __m128i clprod15 = _mm_clmulepi64_si128(temp, key8,
//							0x00);
//					acc = _mm_clmulepi64_si128(acc, key8, 0x10);
//					 __m128i t1 = _mm_xor_si128(x1, clprod1);
//					 __m128i t2 = _mm_xor_si128(clprod2, clprod3);
//					 __m128i t3 = _mm_xor_si128(clprod4, clprod5);
//					 __m128i t4 = _mm_xor_si128(clprod6, clprod7);
//					 __m128i t5 = _mm_xor_si128(clprod8, clprod9);
//					 __m128i t6 = _mm_xor_si128(clprod10, clprod11);
//					 __m128i t7 = _mm_xor_si128(clprod12, clprod13);
//					 __m128i t8 = _mm_xor_si128(clprod14, clprod15);
//
//					 __m128i z1 = _mm_xor_si128(t1, t2);
//					 __m128i z2 = _mm_xor_si128(t3, t4);
//					 __m128i z3 = _mm_xor_si128(t5, t6);
//					 __m128i z4 = _mm_xor_si128(t7, t8);
//					 __m128i Z1 = _mm_xor_si128(z1, z2);
//					 __m128i Z2 = _mm_xor_si128(z3, z4);
//					acc = precompReduction64_si128(
//							_mm_xor_si128(acc, _mm_xor_si128(Z1, Z2)));
//
//				}
//
//			}
//			for (; string + 3 < endstring; string += 4) {
//				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
//				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
//				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
//				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
//				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
//				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
//				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
//				acc = precompReduction64_si128(
//						_mm_xor_si128(acc,
//								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
//										_mm_xor_si128(clprod2, clprod3))));
//			}
//		}
//		for (; string + 1 < endstring; string += 2) {
//			__m128i temp = _mm_lddqu_si128((__m128i *) string);
//			acc = _mm_clmulepi64_si128(acc, key, 0x10);
//			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
//			acc = _mm_xor_si128(clprod1, acc);
//			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
//			acc = precompReduction64_si128(acc);
//		}
//	}
//	if (string < endstring) {
//		const __m128i temp = _mm_set_epi64x(0, *string);
//		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
//		acc = precompReduction64_si128(multi);
//		acc = _mm_xor_si128(acc, temp);
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
//// experimental, expected to be faster than fasthashGaloisFieldPoly64_8
//// but is not
//uint64_t halfhashGaloisFieldPoly64_8(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const uint64_t * randomsource = (const uint64_t *) rs;
//	assert(*randomsource != 0); //otherwise silly
//	const uint64_t * const endstring = string + length;
//	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	if (string + 1 < endstring) {
//		// we start by precomputing the powers of the key
//		__m128i tkey2 = precompReduction64_si128(
//				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
//		// powers of the keys are packed into two registers
//		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
//				_mm_slli_si128(tkey2, 8));
//		if (string + 3 < endstring) {
//			__m128i tkey3 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
//			__m128i tkey4 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
//			__m128i key2 = _mm_xor_si128(
//					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
//					_mm_slli_si128(tkey4, 8));
//			if (string + 7 < endstring) {
//				__m128i tkey5 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey2, tkey3, 0x00));
//				__m128i tkey6 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey3, tkey3, 0x00));
//				__m128i tkey7 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey3, tkey4, 0x00));
//				__m128i tkey8 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey4, tkey4, 0x00));
//				__m128i key3 = _mm_xor_si128(
//						_mm_and_si128(tkey5, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey6, 8));
//				__m128i key4 = _mm_xor_si128(
//						_mm_and_si128(tkey7, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey8, 8));
//
//				/**
//				 * For the half multiplication thing, we want to have the keys organized
//				 * as follows:
//				 * (0,k),(k^2,k^3)...
//				 */
//				// TODO: these keys are probably wrong
//				__m128i hkey1 = _mm_slli_si128(key,8);
//				__m128i hkey2 = _mm_xor_si128(_mm_srli_si128(key,8), _mm_slli_si128(key2,8));
//				__m128i hkey3 = _mm_xor_si128(_mm_srli_si128(key2,8), _mm_slli_si128(key3,8));
//				__m128i hkey4 = _mm_xor_si128(_mm_srli_si128(key3,8), _mm_slli_si128(key4,8));
//				for (; string + 7 < endstring; string += 8) {
//					__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
//					__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
//					__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4)); //a5 a6
//					__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6)); //a7 a8
//					__m128i t1 = _mm_xor_si128(hkey1,temp);
//					__m128i t2 = _mm_xor_si128(hkey2,temp2);
//					__m128i t3 = _mm_xor_si128(hkey3,temp3);
//					__m128i t4 = _mm_xor_si128(hkey4,temp4);
//					__m128i clprod1 = _mm_clmulepi64_si128(t1,t1,0x10);
//					__m128i clprod2 = _mm_clmulepi64_si128(t2,t2,0x10);
//					__m128i clprod3 = _mm_clmulepi64_si128(t3,t3,0x10);
//					__m128i clprod4 = _mm_clmulepi64_si128(t4,t4,0x10);
//					__m128i tacc = _mm_clmulepi64_si128(acc, key4, 0x10); //k^8
//					const __m128i b1 = _mm_xor_si128(clprod2, clprod1);
//					const __m128i b2 = _mm_xor_si128(clprod3, clprod4);
//					const __m128i z1 = _mm_xor_si128(b1, b2);
//					acc = precompReduction64_si128
//							(
//							_mm_xor_si128(tacc, z1));
//				}
//
//			}
//			for (; string + 3 < endstring; string += 4) {
//				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
//				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
//				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
//				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
//				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
//				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
//				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
//				acc = precompReduction64_si128
//						(
//						_mm_xor_si128(acc,
//								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
//										_mm_xor_si128(clprod2, clprod3))));
//			}
//		}
//		for (; string + 1 < endstring; string += 2) {
//			__m128i temp = _mm_lddqu_si128((__m128i *) string);
//			acc = _mm_clmulepi64_si128(acc, key, 0x10);
//			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
//			acc = _mm_xor_si128(clprod1, acc);
//			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
//			acc = precompReduction64_si128(acc);
//		}
//	}
//	if (string < endstring) {
//		const __m128i temp = _mm_set_epi64x(0, *string);
//		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
//		acc = precompReduction64_si128(multi);
//		acc = _mm_xor_si128(acc, temp);
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
//
//uint64_t halfhashGaloisFieldPoly64_16(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const uint64_t * randomsource = (const uint64_t *) rs;
//	assert(*randomsource != 0); //otherwise silly
//	const uint64_t * const endstring = string + length;
//	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//	if (string + 1 < endstring) {
//		// we start by precomputing the powers of the key
//		__m128i tkey2 = precompReduction64_si128(
//				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
//		// powers of the keys are packed into two registers
//		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
//				_mm_slli_si128(tkey2, 8));
//		if (string + 3 < endstring) {
//			__m128i tkey3 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
//			__m128i tkey4 = precompReduction64_si128(
//					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
//			__m128i key2 = _mm_xor_si128(
//					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
//					_mm_slli_si128(tkey4, 8));
//			if (string + 31 < endstring) {
//				__m128i tkey5 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey2, tkey3, 0x00));
//				__m128i tkey6 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey3, tkey3, 0x00));
//				__m128i tkey7 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey3, tkey4, 0x00));
//				__m128i tkey8 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey4, tkey4, 0x00));
//				__m128i tkey9 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey5, tkey4, 0x00));
//				__m128i tkey10 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey5, tkey5, 0x00));
//				__m128i tkey11 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey6, tkey5, 0x00));
//				__m128i tkey12 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey6, tkey6, 0x00));
//				__m128i tkey13 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey6, tkey7, 0x00));
//				__m128i tkey14 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey7, tkey7, 0x00));
//				__m128i tkey15 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey8, tkey7, 0x00));
//				__m128i tkey16 = precompReduction64_si128(
//						_mm_clmulepi64_si128(tkey8, tkey8, 0x00));
//				__m128i key3 = _mm_xor_si128(
//						_mm_and_si128(tkey5, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey6, 8));
//				__m128i key4 = _mm_xor_si128(
//						_mm_and_si128(tkey7, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey8, 8));
//				__m128i key5 = _mm_xor_si128(
//						_mm_and_si128(tkey9, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey10, 8));
//				__m128i key6 = _mm_xor_si128(
//						_mm_and_si128(tkey11, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey12, 8));
//				__m128i key7 = _mm_xor_si128(
//						_mm_and_si128(tkey13, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey14, 8));
//				__m128i key8 = _mm_xor_si128(
//						_mm_and_si128(tkey15, _mm_set_epi64x(0, -1)),
//						_mm_slli_si128(tkey16, 8));
//
//				/**
//				 * For the half multiplication thing, we want to have the keys organized
//				 * as follows:
//				 * (0,k),(k^2,k^3)...
//				 */
//				// TODO: these keys are probably wrong
//				__m128i hkey1 = _mm_slli_si128(key,8);
//				__m128i hkey2 = _mm_xor_si128(_mm_srli_si128(key,8), _mm_slli_si128(key2,8));
//				__m128i hkey3 = _mm_xor_si128(_mm_srli_si128(key2,8), _mm_slli_si128(key3,8));
//				__m128i hkey4 = _mm_xor_si128(_mm_srli_si128(key3,8), _mm_slli_si128(key4,8));
//				__m128i hkey5 = _mm_xor_si128(_mm_srli_si128(key4,8), _mm_slli_si128(key5,8));
//				__m128i hkey6 = _mm_xor_si128(_mm_srli_si128(key5,8), _mm_slli_si128(key6,8));
//				__m128i hkey7 = _mm_xor_si128(_mm_srli_si128(key6,8), _mm_slli_si128(key7,8));
//				__m128i hkey8 = _mm_xor_si128(_mm_srli_si128(key7,8), _mm_slli_si128(key8,8));
//
//				for (; string + 15 < endstring; string += 16) {
//					__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
//					__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
//					__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4)); //a5 a6
//					__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6)); //a7 a8
//					__m128i temp5 = _mm_lddqu_si128((__m128i *) (string + 8)); //a7 a8
//					__m128i temp6 = _mm_lddqu_si128((__m128i *) (string + 10)); //a7 a8
//					__m128i temp7 = _mm_lddqu_si128((__m128i *) (string + 12)); //a7 a8
//					__m128i temp8 = _mm_lddqu_si128((__m128i *) (string + 14)); //a7 a8
//					__m128i t1 = _mm_xor_si128(hkey1,temp);
//					__m128i t2 = _mm_xor_si128(hkey2,temp2);
//					__m128i t3 = _mm_xor_si128(hkey3,temp3);
//					__m128i t4 = _mm_xor_si128(hkey4,temp4);
//					__m128i t5 = _mm_xor_si128(hkey5,temp5);
//					__m128i t6 = _mm_xor_si128(hkey6,temp6);
//					__m128i t7 = _mm_xor_si128(hkey7,temp7);
//					__m128i t8 = _mm_xor_si128(hkey8,temp8);
//					__m128i clprod1 = _mm_clmulepi64_si128(t1,t1,0x10);
//					__m128i clprod2 = _mm_clmulepi64_si128(t2,t2,0x10);
//					__m128i clprod3 = _mm_clmulepi64_si128(t3,t3,0x10);
//					__m128i clprod4 = _mm_clmulepi64_si128(t4,t4,0x10);
//					__m128i clprod5 = _mm_clmulepi64_si128(t5,t5,0x10);
//					__m128i clprod6 = _mm_clmulepi64_si128(t6,t6,0x10);
//					__m128i clprod7 = _mm_clmulepi64_si128(t7,t7,0x10);
//					__m128i clprod8 = _mm_clmulepi64_si128(t8,t8,0x10);
//					__m128i tacc = _mm_clmulepi64_si128(acc, key8, 0x10);
//					const __m128i b1 = _mm_xor_si128(clprod2, clprod1);
//					const __m128i b2 = _mm_xor_si128(clprod3, clprod4);
//					const __m128i b3 = _mm_xor_si128(clprod5, clprod6);
//					const __m128i b4 = _mm_xor_si128(clprod7, clprod8);
//					const __m128i z1 = _mm_xor_si128(b1, b2);
//					const __m128i z2 = _mm_xor_si128(b3, b4);
//					const __m128i Z1 = _mm_xor_si128(z1, z2);
//					acc = precompReduction64_si128
//							(
//							_mm_xor_si128(tacc, Z1));
//				}
//
//			}
//			for (; string + 3 < endstring; string += 4) {
//				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
//				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
//				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
//				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
//				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
//				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
//				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
//				acc = precompReduction64_si128
//						(
//						_mm_xor_si128(acc,
//								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
//										_mm_xor_si128(clprod2, clprod3))));
//			}
//		}
//		for (; string + 1 < endstring; string += 2) {
//			__m128i temp = _mm_lddqu_si128((__m128i *) string);
//			acc = _mm_clmulepi64_si128(acc, key, 0x10);
//			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
//			acc = _mm_xor_si128(clprod1, acc);
//			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
//			acc = precompReduction64_si128(acc);
//		}
//	}
//	if (string < endstring) {
//		const __m128i temp = _mm_set_epi64x(0, *string);
//		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
//		acc = precompReduction64_si128(multi);
//		acc = _mm_xor_si128(acc, temp);
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
//
//uint64_t clmulcacheline(const void* rs, const uint64_t * string,
//		const size_t length) {
//	const __m128i * randomsource = (const __m128i *) rs;
//	// 4 128-bit is a cache line!!!
//	__m128i key1 = _mm_lddqu_si128( randomsource);
//	__m128i key2 = _mm_lddqu_si128( randomsource + 1 );
//	__m128i key3 = _mm_lddqu_si128( randomsource + 2 );
//	__m128i key4 = _mm_lddqu_si128( randomsource + 3 );
//
//	const uint64_t * const endstring = string + length;
//	__m128i acc = _mm_set_epi64x(0, *string);
//	++string;
//
//	for(;string + 7 < endstring;string+=8) {
//		__m128i p1 = _mm_clmulepi64_si128(acc, key1, 0x00);
//
//		__m128i temp1 = _mm_lddqu_si128((__m128i *) string);
//		__m128i p2 = _mm_clmulepi64_si128(temp1, key1, 0x01);
//		__m128i p3 = _mm_clmulepi64_si128(temp1, key2, 0x00);
//
//		__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2));
//		__m128i p4 = _mm_clmulepi64_si128(temp2, key2, 0x01);
//		__m128i p5 = _mm_clmulepi64_si128(temp2, key3, 0x00);
//
//		__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4));
//		__m128i p6 = _mm_clmulepi64_si128(temp3, key3, 0x01);
//		__m128i p7 = _mm_clmulepi64_si128(temp3, key4, 0x00);
//
//		__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6));
//		__m128i p8 = _mm_clmulepi64_si128(temp4, key4, 0x01);
//		__m128i p9 = _mm_slli_si128(temp4, 8);
//
//
//		__m128i q1 = _mm_xor_si128(p1,p2);
//		__m128i q2 = _mm_xor_si128(p3,p4);
//		__m128i q3 = _mm_xor_si128(p5,p6);
//		__m128i q4 = _mm_xor_si128(p7,p8);
//
//		__m128i r1 = _mm_xor_si128(q1,q2);
//		__m128i r2 = _mm_xor_si128(q3,q4);
//		__m128i s1 = _mm_xor_si128(r1,r2);
//		__m128i s2 = _mm_xor_si128(s1,p9);
//		acc = precompReduction64_si128(s2);
//	}
//
//	for (;string + 1 < endstring;string+=2) {
//		__m128i p1 = _mm_clmulepi64_si128(acc, key1, 0x00);
//		__m128i temp1 = _mm_lddqu_si128((__m128i *) string);
//		__m128i p2 = _mm_clmulepi64_si128(temp1, key1, 0x01);
//		__m128i p3 = _mm_slli_si128(temp1, 8);
//		acc = precompReduction64_si128(_mm_xor_si128(_mm_xor_si128(p1,p2),p3));
//	}
//	if(string<endstring){
//		__m128i p1 = _mm_clmulepi64_si128(acc, key1, 0x00);
//		acc = precompReduction64_si128(_mm_xor_si128(p1,_mm_set_epi64x(0, *string)));
//	}
//	return _mm_cvtsi128_si64(acc);
//}
//
