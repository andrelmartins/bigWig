BWPATH:=$(shell Rscript bigWigpath.R)

all: bwrange

bwrange: bwrange.c $(BWPATH)/lib/bigwiglib.a
	gcc -Wall -g -o $@ -I $(BWPATH)/include bwrange.c $(BWPATH)/lib/bigwiglib.a -lm -lz -pthread `pkg-config openssl --libs`

clean:
	rm -f bwrange
