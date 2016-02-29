.PHONY: all clean smhasherpackage-target

all: $(patsubst %.c,%.o,$(wildcard *.c)) smhasherpackage-target

smhasherpackage-target: $(patsubst %.c,%.o,$(wildcard *.c))
	$(MAKE) -C smhasherpackage

%.o: %.c $(wildcard ../../include/*.h)
	$(CC) $(CFLAGS) -o $@ -c $< -I../../../include

clean:
	rm -f *.o *.exe
