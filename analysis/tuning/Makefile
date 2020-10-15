.phony: all clean

%.exe: %.cc ../../include/treehash/*
	$(CXX) $(CXXFLAGS) -I../../include -o $@ $<

clean:
	rm -f *.exe
