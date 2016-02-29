.PHONY: all clean

all: $(patsubst %.c,%.exe,$(wildcard *.c)) $(patsubst %.cc,%.exe,$(wildcard *.cc))

%.exe: %.c $(wildcard ../../include/*.h)
	$(CC) $(CFLAGS) -o $@ $< ../../include/*/*.o -I../../include

%.exe: %.cc $(wildcard ../../include/*.h) $(wildcard ../../include/*.hh)
	$(CXX) $(CXXFLAGS) -o $@ $< ../../include/*/*.o -I../../include

clean:
	rm -f *.o *.exe
