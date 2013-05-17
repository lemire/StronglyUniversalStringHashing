.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

CXXFLAGS =  -O3 -Wall -march=native -Wextra -Wcast-align  

all:  multilinearhashing

multilinearhashing: multilinearhashing.c  
	$(CXX) $(CXXFLAGS) -o multilinearhashing multilinearhashing.c  



clean: 
	rm -f multilinearhashing
