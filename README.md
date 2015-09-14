StronglyUniversalStringHashing
==============================
[![Build Status](https://travis-ci.org/lemire/StronglyUniversalStringHashing.png)](https://travis-ci.org/lemire/StronglyUniversalStringHashing)

Benchmark showing that we can randomly hash strings very quickly with good universality 

This code assumes that you understand random hashing and what it entails.



Further reading:

 Daniel Lemire and Owen Kaser, Faster 64-bit universal hashing using carry-less multiplications, Journal of Cryptographic Engineering (to appear)
 http://arxiv.org/abs/1503.03465 
 
 Owen Kaser and Daniel Lemire, Strongly universal string hashing is fast, Computer Journal (2014) 57 (11): 1624-1638.
 http://arxiv.org/abs/1202.4961




Acknowledgements
==================


Thanks to Nathan Kurz for noticing that GCC 4.7 requires no-tree-vectorize to produce correct results.



Usage
======

    make
    ./benchmark
    ./variablelengthbenchmark

If you plan to use clmul instructions, please run the corresponding
tests:

    make clmulunit
    ./clmulunit


Licensing
==========

In a subdirectory, we have included a modified version of smhasher which is covered under
the MIT license.

In another subdirectory, we have included an implementation of VHASH. It has been put in the
public domain by its authors.

In yet another directory, we have included a C port of CityHash, published under a MIT license
by Google.


Related projects
=================

For a C++ project with similar goals, see:

https://github.com/lemire/fasthashing


