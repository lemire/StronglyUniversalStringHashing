.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

CFLAGS = -std=gnu99  -O3 -Wall -march=native -Wextra -Wcast-align  

all:  multilinearhashing

multilinearhashing: multilinearhashing.c  
	$(CC) $(CFLAGS) -o multilinearhashing multilinearhashing.c  



clean: 
	rm -f multilinearhashing
