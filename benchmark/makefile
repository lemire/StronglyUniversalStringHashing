.phony: all clean

all: $(patsubst %.c,%.exe,$(wildcard *.c)) $(patsubst %.cc,%.exe,$(wildcard *.cc))

%.exe: %.cc $(shell find ../include)
	$(CXX) $(CXXFLAGS) -o $@ $< ../include/*/*.o -I../include  

%.exe: %.c $(shell find ../include)
	$(CC) $(CFLAGS) -o $@ $< ../include/*/*.o -I../include  -lstdc++

clean:
	rm -f *.o *.exe
