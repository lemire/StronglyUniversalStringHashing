.PHONY: all clean

all: $(patsubst %.c,%.o,$(wildcard *.c))

%.o: %.c *.h
	$(CC) $(CFLAGS) -w -c $<

clean:
	rm -f *.o
