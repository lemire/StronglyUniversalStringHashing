.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

CFLAGS = -std=gnu99 -funroll-loops -O3 -mavx  -march=native   

all: clmulunit  benchmark City.o

benchmark: src/benchmark.c include/*.h City.o 
	$(CC) $(CFLAGS) -o benchmark src/benchmark.c City.o -Iinclude -ICity

variablelengthbenchmark: src/variablelengthbenchmark.c include/*.h City.o 
	$(CC) $(CFLAGS) -o variablelenthbenchmark src/variablelengthbenchmark.c City.o -Iinclude -ICity


clmulunit: src/clmulunit.c include/*.h
	$(CC) $(CFLAGS) -o clmulunit src/clmulunit.c -Iinclude 

City.o: City/City.c City/City.h
	$(CC) $(CFLAGS) -c City/City.c -ICity
 
clean: 
	rm -f multilinearhashing clmulunit *.o
