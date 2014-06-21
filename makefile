.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

CFLAGS = -std=gnu99 -funroll-loops -O3 -mavx  -march=native   

all: clmulunit  benchmark

benchmark: src/benchmark.c include/*.h 
	$(CC) $(CFLAGS) -o benchmark src/benchmark.c -Iinclude 

clmulunit: src/clmulunit.c include/*.h
	$(CC) $(CFLAGS) -o clmulunit src/clmulunit.c -Iinclude 


clean: 
	rm -f multilinearhashing clmulunit *.o
