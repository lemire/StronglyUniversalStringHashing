.phony: all City-target SipHash-target VHASH-target PMP-target clean

CFLAGS = $(FLAGS) -fPIC
CDEBUGFLAGS = $(DEBUGFLAGS) -fPIC
CXXFLAGS = $(FLAGS) -fPIC
CXXDEBUGFLAGS = $(DEBUGFLAGS) -fPIC
export

all: City-target SipHash-target VHASH-target PMP-target umash/umash.o

%-target:
	$(MAKE) -C $*

clean:
	$(MAKE) -C City clean
	$(MAKE) -C SipHash clean
	$(MAKE) -C VHASH clean
	$(MAKE) -C PMP clean
	rm -f umash/umash.o
