.phony: all iaca-target tuning-target clean

all: iaca-target tuning-target

%-target:
	$(MAKE) -C $*

clean:
	$(MAKE) -C iaca clean
	$(MAKE) -C tuning clean
