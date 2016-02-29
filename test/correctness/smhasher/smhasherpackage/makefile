.PHONY: all clean

all: $(patsubst %.c,%.o,$(wildcard *.c)) $(patsubst %.cpp,%.o,$(wildcard *.cpp))

%.o: %.c *.h
	$(CC) $(CFLAGS) -w -c $<

%.o: %.cpp *.h
	$(CXX) $(CXXFLAGS) -w -c $<

smhasher.exe: $(patsubst %.c,%.o,$(wildcard *.c)) $(patsubst %.cpp,%.o,$(wildcard *.cpp))
	$(CXX) $(CXXFLAGS) -o $@ $^ ../*.o

clean:
	rm -f *.o
