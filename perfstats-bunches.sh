#!/bin/bash

bunches=( cycles,instructions,cache-references,cache-misses,branches branch-misses,bus-cycles,L1-dcache-loads,L1-dcache-load-misses L1-dcache-store-misses,L1-dcache-prefetch-misses,L1-icache-load-misses dTLB-load-misses,dTLB-store-misses,iTLB-load-misses r0203,r0803,r0105,r0205 r0107,r0108,r0208,r0408,r0e08 r1008,r2008,r4008,r6008,r8008 )

#many more to be added

for b in  ${bunches[*]} ; do
#    echo bunch $b
    ./perfstats.sh $b
    echo
done



