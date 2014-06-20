StronglyUniversalStringHashing
==============================

Benchmark showing the we can randomly hash strings very quickly with good strong universality 

This code assumes that you understand random hashing and what it entails.

 Reference: Owen Kaser and Daniel Lemire, Strongly universal string hashing is fast, Computer Journal (to appear)
 http://arxiv.org/abs/1202.4961


For a C++ equivalent, see:

https://github.com/lemire/fasthashing


Thanks to Nathan Kurz for noticing that GCC 4.7 requires no-tree-vectorize to produce correct results.



Usage
======


If you plan to use clmul instructions, please run the corresponding
tests:

make clmulunit
./clmulunit