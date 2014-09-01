.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

CXXFLAGS = -O2 -mavx -march=native  -Wall 

all:  smhasher

smhasher: *.h *.cpp *.c cl3264.o         vhash4smhasher.o  vmac.o rijndael-alg-fst.o
	$(CXX) $(CXXFLAGS)  -o smhasher *cpp *.c cl3264.o  vhash4smhasher.o vmac.o rijndael-alg-fst.o


clean: 
	rm -f smhasher *.o
