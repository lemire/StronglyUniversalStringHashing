#ifndef PMP_C_WRAPPER
#define PMP_C_WRAPPER


/*
C interface to the PMP library.

TODO: implement seeding and 64-bit hashing
*/
#ifdef __cplusplus
extern "C" {
#endif

/* Hash the string made of length characters, returns a 64-bit value.*/
uint64_t pmp64_hash( const unsigned char* chars, size_t length);



#ifdef __cplusplus
} // extern "C"
#endif

#endif //PMP_C_WRAPPER
