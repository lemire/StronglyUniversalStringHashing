.phony: all clean

all: $(patsubst %.c,%.exe,$(wildcard *.c)) $(patsubst %.cc,%.exe,$(wildcard *.cc))

%.exe: %.cc
	$(CXX) $(CXXFLAGS) -o $@ $< ../include/*/*.o -I../include

%.exe: %.c
	$(CC) $(CFLAGS) -o $@ $< ../include/*/*.o -I../include

clean:
	rm -f *.o *.exe
