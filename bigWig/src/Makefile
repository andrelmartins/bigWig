include bigwig.mk
include jkweb/inc/common.mk

# Unix specific linker flags
L += -lz

# build our MAKEFLAGS
# R needs spaces to be escaped ...
#
MFLGS := "CFLAGS=-I./jkweb/inc ${COPT} ${CFLAGS}"
# define a space
empty:=
space:= $(empty) $(empty)
MFLGSESC := $(subst $(space),\$(space),$(MFLGS))

all: bigwig.so bigwiglib.a bigwiglib.so

jkweb.a:
	[ -d jkweb/lib/${MACHTYPE} ] || mkdir jkweb/lib/${MACHTYPE}
	$(MAKE) -C jkweb/lib
	ln -s jkweb/lib/${MACHTYPE}/jkweb.a jkweb.a	

bigwig.so: prepare jkweb.a bigWig_R.c utils_R.c bwgExtra.h bwgExtra.c bw_query.h bw_query.c bw_base.h bw_base.c bigWig_R2.c
	@MAKEFLAGS=${MFLGSESC} $(R_HOME)/bin/R CMD SHLIB -o $@ bigWig_R.c utils_R.c bwgExtra.c  bw_base.c bw_query.c bigWig_R2.c jkweb.a ${L}

bigwiglib.a: jkweb.a bw_query.h bw_query.c bw_base.h bw_base.c bigwiglib.c
	$(CC) -I./jkweb/inc ${COPT} ${CFLAGS} -c -o bigwiglib_base.o bw_base.c
	$(CC) -I./jkweb/inc ${COPT} ${CFLAGS} -c -o bigwiglib_query.o bw_query.c
	$(CC) -I./jkweb/inc ${COPT} ${CFLAGS} -c -o bigwiglib.o bigwiglib.c
	ar rcs $@ bigwiglib.o bigwiglib_base.o bigwiglib_query.o `find jkweb/ -name "*.o"`
	[ -d ../inst/lib ] || mkdir ../inst/lib
	cp bigwiglib.a ../inst/lib

bigwiglib.so: jkweb.a bw_query.h bw_query.c bw_base.h bw_base.c bigwiglib.c
	$(CC) -shared -I./jkweb/inc ${COPT} ${CFLAGS} -o $@ bw_base.c bw_query.c bigwiglib.c jkweb.a ${L} ${LDFLAGS}

clean:
	rm -f bigwig.so bigwiglib.a bigwiglib.so *.o jkweb.a
	$(MAKE) -C jkweb/lib clean
	rm -f ../inst/lib/bigwiglib.a
	rm -Rf bigwiglib.so.dSYM

tgtclean:
	rm -f bigwig.so bigwiglib.a bigwiglib.so *.o

.PHONY: prepare
