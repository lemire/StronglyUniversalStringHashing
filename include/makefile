.phony: all City-target SipHash-target VHASH-target PMP-target clean

all: City-target SipHash-target VHASH-target PMP-target

%-target:
	$(MAKE) -C $*

clean:
	$(MAKE) -C City clean
	$(MAKE) -C SipHash clean
	$(MAKE) -C VHASH clean
	$(MAKE) -C PMP clean
