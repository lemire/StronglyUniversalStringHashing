.PHONY: all clean

all: $(patsubst %.c,%.o,$(wildcard *.c))

%.o: %.c $(wildcard ./*.h)
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o
