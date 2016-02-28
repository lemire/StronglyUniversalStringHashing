.PHONY: all clean

all: $(patsubst %.c,%.exe,$(wildcard *.c)) $(patsubst %.cc,%.exe,$(wildcard *.cc))

%.exe: %.c $(wildcard ../../include/*.h)
	$(CC) $(CFLAGS) -o $@ $< -I../../include

%.exe: %.cc $(wildcard ../../include/*.h) $(wildcard ../../include/*.hh)
	$(CXX) $(CXXFLAGS) -o $@ $< ../../City.o ../../siphash24.o ../../vmac.o ../../rijndael-alg-fst.o -I../../include -I../../City -I../../SipHash -I../../VHASH

clean:
	rm -f *.o *.exe
