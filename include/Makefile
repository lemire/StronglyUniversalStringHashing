.phony: all City-target SipHash-target VHASH-target PMP-target clean

CFLAGS = $(FLAGS) -fPIC  -std=gnu11
CDEBUGFLAGS = $(DEBUGFLAGS) -fPIC  -std=gnu11
CXXFLAGS = $(FLAGS) -fPIC -std=c++11
CXXDEBUGFLAGS = $(DEBUGFLAGS) -fPIC -std=c++11
export

all: City-target SipHash-target VHASH-target PMP-target umash/README.md umash/umash.o
umash/README.md:
	git submodule update --init --recursive

umash/umash.o: umash/README.md

%-target:
	$(MAKE) -C $*

clean:
	$(MAKE) -C City clean
	$(MAKE) -C SipHash clean
	$(MAKE) -C VHASH clean
	$(MAKE) -C PMP clean
	rm -f umash/umash.o
