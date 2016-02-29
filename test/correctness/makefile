.PHONY: all clean smhasher-target

all: $(patsubst %.c,%.exe,$(wildcard *.c)) \
     $(patsubst %.cc,%.exe,$(wildcard *.cc)) smhasher-target

%.exe: %.c $(wildcard ../../include/*.h)
	$(CC) $(CFLAGS) -o $@ $< ../../include/*/*.o -I../../include  -lstdc++

%.exe: %.cc $(wildcard ../../include/*.h) $(wildcard ../../include/*.hh)
	$(CXX) $(CXXFLAGS) -o $@ $< ../../include/*/*.o -I../../include

smhasher-target:
	$(MAKE) -C smhasher

%.o: %.c $(wildcard ../../include/*.h)
	$(CC) $(CFLAGS) -o $@ -c $< -I../../include

smhasher.exe: smhasherpackage cl3264.o vhash4smhasher.o
	$(MAKE) -C smhasherpackage
	$(CXX) $(FLAGS) -o $@ smhasherpackage/*.o cl3264.o vhash4smhasher.o \
            ../../include/*/*.o

clean:
	rm -f *.o *.exe
	$(MAKE) -C smhasher clean
