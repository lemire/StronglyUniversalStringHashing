.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

CFLAGS = -std=gnu99 -funroll-loops -O3 -mavx  -march=native   

all:  multilinearhashing

multilinearhashing: src/multilinearhashing.c include/*.h 
	$(CC) $(CFLAGS) -o multilinearhashing src/multilinearhashing.c -Iinclude 

clmulunit: src/clmulunit.c include/clmul.h
	$(CC) $(CFLAGS) -o clmulunit src/clmulunit.c -Iinclude 


clean: 
	rm -f multilinearhashing clmulunit *.o
