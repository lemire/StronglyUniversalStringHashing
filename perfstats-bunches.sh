#!/bin/bash

bunches=( cycles,instructions,cache-references,cache-misses,branches branch-misses,bus-cycles,L1-dcache-loads,L1-dcache-load-misses L1-dcache-store-misses,L1-dcache-prefetch-misses,L1-icache-load-misses dTLB-load-misses,dTLB-store-misses,iTLB-load-misses r0203,r0803,r0105,r0205 r0107,r0108,r0208,r0408 r0e08,r1008,r2008,r4008 r6008,r8008,r030d,r010e r2124,r2224,r2424,rf824 r3f24,r0148,r0149,r2049 r6049,r0151,r0458,r0858 r015e,r0163,r0263,r0279 r0280,r0185,r1085,r0187 r0487,rff88,rff89,r019c r01a2,r02a3,r08c1,r10c1 r40c1,r20c3,r01c4,r02c4 r04c4,r08c4,r10c4,r20c4 r40c4,r1eca,r1fe6,r07f1 r04a2,r08a2,r10a2)

#many more to be added

for b in  ${bunches[*]} ; do
#    echo bunch $b
    ./perfstats.sh $b
    echo
done



