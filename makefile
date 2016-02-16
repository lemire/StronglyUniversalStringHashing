.SUFFIXES:

.phony: all clean smhasherpackage-submake

FLAGS = -ggdb -O2 -mavx -mavx2 -march=native -Wall -Wextra -Wstrict-overflow -Wstrict-aliasing -funroll-loops
DEBUGFLAGS = $(FLAGS) -ggdb3 -O0 -fno-unroll-loops -fsanitize=undefined
CFLAGS = $(FLAGS) -std=gnu99
CDEBUGFLAGS = $(DEBUGFLAGS) -std=gnu99
CXXFLAGS = $(FLAGS) -std=c++0x
CXXDEBUGFLAGS = $(DEBUGFLAGS) -std=c++0x
export

all: clmulunit classicvariablelengthbenchmark variablelengthbenchmark variablelengthbenchmark-debug benchmark benchmark64bitreductions uniformsanity smhasher benchmark128bitmultiplication benchmark128bitpolyhashing boosted-treehash-params.exe

nhvsclnh.o: src/nhvsclnh.c
	$(CC) $(CFLAGS)  -c  src/nhvsclnh.c

uniformsanity: src/uniform_sanity.c include/*.h City.o vmac.o siphash24.o
	$(CC) $(CFLAGS)  -o uniformsanity src/uniform_sanity.c City.o siphash24.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -ISipHash -IVHASH 

benchmark: src/benchmark.c include/*.h City.o siphash24.o vmac.o 
	$(CC) $(CFLAGS) -o benchmark src/benchmark.c City.o siphash24.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -ISipHash -IVHASH 

benchmark-debug: src/benchmark.c include/*.h City.o siphash24.o vmac.o 
	$(CC) $(CDEBUGFLAGS) -o benchmark-debug src/benchmark.c City.o siphash24.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -ISipHash -IVHASH 

benchmark128bitmultiplication: src/benchmark128bitmultiplication.c include/clmul.h  
	$(CC) $(CFLAGS) -o benchmark128bitmultiplication src/benchmark128bitmultiplication.c   -Iinclude 

benchmark128bitpolyhashing: src/benchmark128bitpolyhashing.c include/clmul.h  
	$(CC) $(CFLAGS) -o benchmark128bitpolyhashing src/benchmark128bitpolyhashing.c   -Iinclude 

benchmark64bitreductions: src/benchmark64bitreductions.c include/clmul.h  
	$(CC) $(CFLAGS) -o benchmark64bitreductions src/benchmark64bitreductions.c   -Iinclude 

classicvariablelengthbenchmark: src/classicvariablelengthbenchmark.c include/*.h City.o siphash24.o vmac.o 
	$(CC) $(CFLAGS) -o classicvariablelengthbenchmark src/classicvariablelengthbenchmark.c City.o siphash24.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -ISipHash -IVHASH

variablelengthbenchmark: src/variablelengthbenchmark.cc include/*.h include/treehash/*.hh City.o siphash24.o vmac.o 
	$(CXX) $(CXXFLAGS) -o variablelengthbenchmark src/variablelengthbenchmark.cc City.o siphash24.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -ISipHash -IVHASH 

variablelengthbenchmark-debug: src/variablelengthbenchmark.cc include/*.h include/treehash/*.hh City.o siphash24.o vmac.o 
	$(CXX) $(CXXDEBUGFLAGS) -o $@ src/variablelengthbenchmark.cc City.o siphash24.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -ISipHash -IVHASH 

boosted-treehash-params.exe: src/boosted-treehash-params.cc include/*.h include/treehash/*.hh
	$(CXX) $(CXXFLAGS) -Iinclude -o $@ $<

clmulunit: src/clmulunit.c include/*.h
	$(CC) $(CFLAGS) -o clmulunit src/clmulunit.c -Iinclude 

City.o: City/City.c City/City.h
	$(CC) $(CFLAGS) -c City/City.c -ICity 
 
siphash24.o: SipHash/siphash24.c SipHash/siphash24.h
	$(CC) $(CFLAGS) -c SipHash/siphash24.c -ISipHash 


rijndael-alg-fst.o: VHASH/rijndael-alg-fst.c  VHASH/rijndael-alg-fst.h 
	$(CC) $(CFLAGS) -c VHASH/rijndael-alg-fst.c -IVHASH 

cl3264.o:	src/cl3264.c include/*.h
	$(CC) $(CFLAGS) -c src/cl3264.c -Iinclude

vhash4smhasher.o:	src/vhash4smhasher.c include/*.h
	$(CC) $(CFLAGS) -c src/vhash4smhasher.c -Iinclude -IVHASH 

vmac.o: rijndael-alg-fst.o VHASH/vmac.c VHASH/vmac.h
	$(CC) $(CFLAGS) -c VHASH/vmac.c -IVHASH 

smhasher: $(wildcard smhasherpackage/*.h) $(wildcard smhasherpackage/*.c) $(wildcard smhasherpackage/*.cpp) cl3264.o vhash4smhasher.o vmac.o rijndael-alg-fst.o
	$(MAKE) -C smhasherpackage
	$(CXX) $(FLAGS) -o smhasher smhasherpackage/*.o cl3264.o vhash4smhasher.o vmac.o rijndael-alg-fst.o

clean:
	$(MAKE) -C smhasherpackage clean
	rm -f multilinearhashing classicvariablelengthbenchmark variablelengthbenchmark benchmark benchmark64bitreductions clmulunit uniformsanity smhasher variablelenthbenchmark  *.o
