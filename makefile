.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

CFLAGS = -std=gnu99 -ggdb  -O2 -mavx  -march=native   

all: clmulunit variablelengthbenchmark  benchmark benchmark64bitreductions uniformsanity

uniformsanity: src/uniform_sanity.c include/*.h City.o vmac.o
	$(CC) $(CFLAGS)  -o uniformsanity src/uniform_sanity.c City.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -IVHASH 

benchmark: src/benchmark.c include/*.h City.o vmac.o 
	$(CC) $(CFLAGS) -o benchmark src/benchmark.c City.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity -IVHASH 

benchmark64bitreductions: src/benchmark64bitreductions.c include/clmul.h  
	$(CC) $(CFLAGS) -o benchmark64bitreductions src/benchmark64bitreductions.c   -Iinclude 

variablelengthbenchmark: src/variablelengthbenchmark.c include/*.h City.o vmac.o 
	$(CC) $(CFLAGS) -o variablelengthbenchmark src/variablelengthbenchmark.c City.o vmac.o rijndael-alg-fst.o  -Iinclude -ICity  -IVHASH 

clmulunit: src/clmulunit.c include/*.h
	$(CC) $(CFLAGS) -o clmulunit src/clmulunit.c -Iinclude 

City.o: City/City.c City/City.h
	$(CC) $(CFLAGS) -c City/City.c -ICity 
 
rijndael-alg-fst.o: VHASH/rijndael-alg-fst.c  VHASH/rijndael-alg-fst.h 
	$(CC) $(CFLAGS) -c VHASH/rijndael-alg-fst.c -IVHASH 

vmac.o: rijndael-alg-fst.o VHASH/vmac.c VHASH/vmac.h
	$(CC) $(CFLAGS) -c VHASH/vmac.c -IVHASH 
clean: 
	rm -f multilinearhashing variablelengthbenchmark benchmark benchmark64bitreductions clmulunit uniformsanity *.o
