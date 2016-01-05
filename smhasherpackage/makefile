.PHONY: all

all: $(patsubst %.c,%.o,$(wildcard *.c)) $(patsubst %.cpp,%.o,$(wildcard *.cpp))

%.o: %.c *.h
	$(CC) $(CFLAGS) -w -c $<

%.o: %.cpp *.h
	$(CXX) $(CXXFLAGS) -w -c $<

clean:
	rm -f *.o
