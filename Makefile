.SUFFIXES:

.phony: all clean analysis-target test-target benchmark-target

FLAGS = -ggdb -O2 -mavx -mavx2 -march=native -Wall -Wextra -Wstrict-overflow \
        -Wstrict-aliasing -funroll-loops -fno-strict-aliasing
DEBUGFLAGS = $(FLAGS) -ggdb3 -O0 -fno-unroll-loops -fsanitize=undefined
CFLAGS = $(FLAGS) -std=gnu11
CDEBUGFLAGS = $(DEBUGFLAGS) -std=gnu11
CXXFLAGS = $(FLAGS) -std=c++11
CXXDEBUGFLAGS = $(DEBUGFLAGS) -std=c++11
export

all: test-target benchmark-target

analysis-target:
	$(MAKE) -C analysis

test-target:
	$(MAKE) -C include
	$(MAKE) -C test

benchmark-target:
	$(MAKE) -C include
	$(MAKE) -C benchmark

%.exe: src/%.cc $(shell find include -iname '*.h') $(shell find include -iname '*.c')
	$(MAKE) -C include
	$(CXX) $(CXXFLAGS) -o $@ $< include/*/*.o -Iinclude

%.exe: src/%.c $(shell find include -iname '*.h') $(shell find include -iname '*.c')
	$(MAKE) -C include
	$(CC) $(CFLAGS) -o $@ $< include/*/*.o -Iinclude

%.o: src/%.c $(shell find include -iname '*.h') $(shell find include -iname '*.c')
	$(MAKE) -C include
	$(CC) $(CFLAGS) -o $@ -c $< -Iinclude

variablelengthbenchmark-unaligned.exe: variablelengthbenchmark.exe
	ln -sf variablelengthbenchmark.exe variablelengthbenchmark-unaligned.exe

clean:
	$(MAKE) -C analysis clean
	$(MAKE) -C test clean
	$(MAKE) -C benchmark clean
	$(MAKE) -C include clean
