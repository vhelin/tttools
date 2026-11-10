.PHONY: all test clean

all:
	$(MAKE) -C deflateTT
	$(MAKE) -C inflateTT
	$(MAKE) -C inflateTT-MP
	$(MAKE) -C tests tests.exe

test: all
	$(MAKE) -C tests

clean:
	$(MAKE) -C deflateTT clean
	$(MAKE) -C inflateTT clean
	$(MAKE) -C inflateTT-MP clean
	$(MAKE) -C tests clean
