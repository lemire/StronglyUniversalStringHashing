StronglyUniversalStringHashing
==============================
[![Build Status](https://travis-ci.org/lemire/StronglyUniversalStringHashing.png)](https://travis-ci.org/lemire/StronglyUniversalStringHashing)

Very fast [universal hash
families](https://en.wikipedia.org/wiki/Universal_hashing) on strings.

Sample results on a regular x64 (Skylake) processor:
```
Google's City                        CPU cycle/byte = 0.216047 	 
64-bit VHASH                         CPU cycle/byte = 0.215097 	 
64-bit CLHASH                        CPU cycle/byte = 0.091786 	
SipHash                              CPU cycle/byte = 1.414069
```


This software includes the experimental code from the following research papers... 

* Daniel Lemire and Owen
Kaser, [Faster 64-bit universal hashing using carry-less
multiplications](http://arxiv.org/abs/1503.03465), Journal of Cryptographic Engineering 6 (3), 2016.
* Owen Kaser and Daniel Lemire, [Strongly universal string hashing is fast](http://arxiv.org/abs/1202.4961), Computer Journal (2014) 57 (11): 1624-1638.

Acknowledgements
==================

Thanks to Nathan Kurz for noticing that GCC 4.7 requires
no-tree-vectorize to produce correct results.

Usage
======

To test speed:

    make benchmark-target
    # disable some processor features that add noise to benchmarks:
    cd scripts/; sudo ./master.sh; cd ..
    ./benchmark/benchmark.exe
    ./benchmark/variablelengthbenchmark.exe

To test correctness of hash functions using PCLMULQDQ:

    make test-target
    for test in ./test/correctness/*.exe; do $test; done

Or more simply...

    ./run_unit.sh

Related projects
=================

There is a very simple clhash library in C: https://github.com/lemire/clhash

For a project with similar goals, see: https://github.com/lemire/fasthashing

Licensing
==========

This project is licenced as described in the LICENSE file, with the
following exceptions for code written by other authors:

  * smhasher and CityHash are MIT licensed.

  * VHASH and siphash are public domain.


