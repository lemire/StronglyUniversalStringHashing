StronglyUniversalStringHashing
==============================
[![Build Status](https://travis-ci.org/lemire/StronglyUniversalStringHashing.png)](https://travis-ci.org/lemire/StronglyUniversalStringHashing)

Very fast [universal hash
families](https://en.wikipedia.org/wiki/Universal_hashing) on strings.

This code includes the experimental code from [Daniel Lemire and Owen
Kaser, "Faster 64-bit universal hashing using carry-less
multiplications", Journal of Cryptographic Engineering (to
appear)](http://arxiv.org/abs/1503.03465) and [Owen Kaser and Daniel
Lemire, "Strongly universal string hashing is fast", Computer Journal
(2014) 57 (11): 1624-1638](http://arxiv.org/abs/1202.4961).

Acknowledgements
==================

Thanks to Nathan Kurz for noticing that GCC 4.7 requires
no-tree-vectorize to produce correct results.

Usage
======

To test speed:

    make
    ./benchmark
    ./variablelengthbenchmark

To test correctness of hash functions using PCLMULQDQ:

    make clmulunit hashunit
    ./clmulunit ; ./hashunit

Licensing
==========

In a subdirectory, we have included a modified version of smhasher
which is covered under the MIT license.

In another subdirectory, we have included an implementation of
VHASH. It has been put in the public domain by its authors.

In yet another directory, we have included a C port of CityHash,
published under a MIT license by Google.

Related projects
=================

For a project with similar goals, see:

https://github.com/lemire/fasthashing