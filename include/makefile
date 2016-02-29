.phony: all City-target SipHash-target VHASH-target clean

all: City-target SipHash-target VHASH-target

%-target:
	$(MAKE) -C $*

clean:
	$(MAKE) -C City clean
	$(MAKE) -C SipHash clean
	$(MAKE) -C VHASH clean
