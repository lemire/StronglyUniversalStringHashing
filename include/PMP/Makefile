.PHONY: all clean

all: $(patsubst %.cpp,%.o,$(wildcard *.cpp))

%.o: %.cpp $(wildcard ./*.h)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o
