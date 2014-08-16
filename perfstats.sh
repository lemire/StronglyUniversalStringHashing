#!/bin/bash

#usage: perfstats.sh  stat1,stat2,stat3,....   

perf stat -e \{$1\} ./variablelengthbenchmark 1 512 512 2> junk1 >/dev/null
egrep "^\s+[0-9]" junk1 | sed -r "s/^ +([^ ]+) ([^ ]+) .*/\2 \1/"  > junk1a

perf stat -e \{$1\} ./variablelengthbenchmark 2 512 512 2> junk2 >/dev/null
egrep "^\s+[0-9]" junk2 | sed -r "s/^ +([^ ]+) ([^ ]+) .*/\1/"  > junk2a

perf stat -e \{$1\} ./variablelengthbenchmark 4 512 512 2> junk4 >/dev/null
egrep "^\s+[0-9]" junk4 | sed -r "s/^ +([^ ]+) ([^ ]+) .*/\1/"  > junk4a

perf stat -e \{$1\} ./variablelengthbenchmark 8 512 512 2> junk8 >/dev/null
egrep "^\s+[0-9]" junk8 | sed -r "s/^ +([^ ]+) ([^ ]+) .*/\1/"  > junk8a

echo stat VHASH PCMUL City datagen
paste junk1a junk2a junk4a junk8a



