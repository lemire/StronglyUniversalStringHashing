
#include "PMP_Multilinear_64.h"
#include "PMP_C_wrapper.h"

static PMP_Multilinear_Hasher_64 pmp64;


#ifdef __cplusplus
extern "C" {
#endif



/* Hash the string made of length characters, returns a 64-bit value.*/
uint64_t pmp64_hash( const unsigned char* chars, size_t length) {
    return pmp64.hash(chars,length);
}


#ifdef __cplusplus
} // extern "C"
#endif
