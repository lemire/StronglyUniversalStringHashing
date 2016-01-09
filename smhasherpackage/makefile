# makefile will be updated when make did something. This serves as a
# proxy to the upstream makefile that something in this directory has
# changed.
makefile: $(patsubst %.c,%.o,$(wildcard *.c)) $(patsubst %.cpp,%.o,$(wildcard *.cpp))
	touch makefile

%.o: %.c *.h
	$(CC) $(CFLAGS) -w -c $<

%.o: %.cpp *.h
	$(CXX) $(CXXFLAGS) -w -c $<

clean:
	rm -f *.o
