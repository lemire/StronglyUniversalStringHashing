.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

CFLAGS = -std=gnu99 -funroll-loops -O3  -mavx -march=native   

all:  multilinearhashing

multilinearhashing: multilinearhashing.c clmul.h 
	$(CC) $(CFLAGS) -o multilinearhashing multilinearhashing.c  

multilinearhashingcl: multilinearhashing.c  
	$(CC) $(CFLAGS) -mpclmul -o multilinearhashingcl multilinearhashing.c  

clean: 
	rm -f multilinearhashing
