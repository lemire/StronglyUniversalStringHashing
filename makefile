.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

CFLAGS = -std=gnu99 -ggdb  -O2 -mavx  -march=native -Wall -Wextra -funroll-loops
CXXFLAGS = -O2 -mavx  -march=native -std=c++0x

all: clmulunit variablelengthbenchmark  benchmark benchmark64bitreductions uniformsanity smhasher benchmark128bitmultiplication benchmark128bitpolyhashing

nhvsclnh.o: src/nhvsclnh.c
	$(CC) $(CFLAGS)  -c  src/nhvsclnh.c

uniformsanity: src/uniform_sanity.c include/*.h City.o vmac.o siphash24.o
	$(CC) $(CFLAGS)  -o uniformsanity src/uniform_sanity.c City.o siphash24.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -ISipHash -IVHASH 

benchmark: src/benchmark.c include/*.h City.o siphash24.o vmac.o 
	$(CC) $(CFLAGS) -o benchmark src/benchmark.c City.o siphash24.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -ISipHash -IVHASH 

benchmark128bitmultiplication: src/benchmark128bitmultiplication.c include/clmul.h  
	$(CC) $(CFLAGS) -o benchmark128bitmultiplication src/benchmark128bitmultiplication.c   -Iinclude 

benchmark128bitpolyhashing: src/benchmark128bitpolyhashing.c include/clmul.h  
	$(CC) $(CFLAGS) -o benchmark128bitpolyhashing src/benchmark128bitpolyhashing.c   -Iinclude 

benchmark64bitreductions: src/benchmark64bitreductions.c include/clmul.h  
	$(CC) $(CFLAGS) -o benchmark64bitreductions src/benchmark64bitreductions.c   -Iinclude 

variablelengthbenchmark: src/variablelengthbenchmark.cc include/*.h City.o siphash24.o vmac.o 
	$(CXX) $(CXXFLAGS) -o variablelengthbenchmark src/variablelengthbenchmark.cc City.o siphash24.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -ISipHash -IVHASH 

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

smhasher: smhasherpackage/*.h smhasherpackage/*.cpp smhasherpackage/*.c cl3264.o         vhash4smhasher.o  vmac.o rijndael-alg-fst.o
	$(CXX) $(CXXFLAGS)  -o smhasher -x c++ smhasherpackage/*.cpp -x c smhasherpackage/*.c -x none cl3264.o  vhash4smhasher.o vmac.o rijndael-alg-fst.o -Ismhasherpackage

clean: 
	rm -f multilinearhashing variablelengthbenchmark benchmark benchmark64bitreductions clmulunit uniformsanity smhasher variablelenthbenchmark  *.o
