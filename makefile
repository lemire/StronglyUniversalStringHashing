.SUFFIXES:

.phony: all clean analysis-target test-target

FLAGS = -ggdb -O2 -mavx -mavx2 -march=native -Wall -Wextra -Wstrict-overflow \
        -Wstrict-aliasing -funroll-loops
DEBUGFLAGS = $(FLAGS) -ggdb3 -O0 -fno-unroll-loops -fsanitize=undefined
CFLAGS = $(FLAGS) -std=gnu99
CDEBUGFLAGS = $(DEBUGFLAGS) -std=gnu99
CXXFLAGS = $(FLAGS) -std=c++0x
CXXDEBUGFLAGS = $(DEBUGFLAGS) -std=c++0x
export

all: test-target variablelengthbenchmark.exe benchmark.exe \
     benchmark64bitreductions.exe benchmark128bitmultiplication.exe \
     benchmark128bitpolyhashing.exe smhasher boosted-treehash-params.exe

analysis-target:
	$(MAKE) -C analysis

test-target:
	$(MAKE) -C include
	$(MAKE) -C test

%.exe: src/%.cc $(shell find include -iname '*.h') $(shell find include -iname '*.c')
	$(MAKE) -C include
	$(CXX) $(CXXFLAGS) -o $@ $< include/*/*.o -Iinclude

%.exe: src/%.c $(shell find include -iname '*.h') $(shell find include -iname '*.c')
	$(MAKE) -C include
	$(CC) $(CFLAGS) -o $@ $< include/*/*.o -Iinclude

%.o: src/%.c $(shell find include -iname '*.h') $(shell find include -iname '*.c')
	$(MAKE) -C include
	$(CC) $(CFLAGS) -o $@ -c $< -Iinclude

variablelengthbenchmark-unaligned: variablelengthbenchmark
	ln -sf variablelengthbenchmark variablelengthbenchmark-unaligned

boosted-treehash-params.exe: src/boosted-treehash-params.cc include/*.h \
                             include/treehash/*.hh
	$(CXX) $(CXXFLAGS) -Iinclude -o $@ $<

smhasher: $(wildcard smhasherpackage/*.h) $(wildcard smhasherpackage/*.c) \
          $(wildcard smhasherpackage/*.cpp) cl3264.o vhash4smhasher.o
	$(MAKE) -C include
	$(MAKE) -C smhasherpackage
	$(CXX) $(FLAGS) -o smhasher smhasherpackage/*.o cl3264.o \
            vhash4smhasher.o include/*/*.o

clean:
	$(MAKE) -C smhasherpackage clean
	$(MAKE) -C analysis clean
	$(MAKE) -C test clean
	rm -f multilinearhashing classicvariablelengthbenchmark \
            variablelengthbenchmark benchmark benchmark64bitreductions \
            smhasher variablelenthbenchmark *.o
