include jkweb/inc/common.mk

ifeq (${MACHTYPE},)
	export MACHTYPE=$(shell echo "cat(R.version[['arch']], '\\\\n', sep='')" | R -q --no-save | grep -v ">")
endif

CC:=$(shell R CMD config CC)
LDFLAGS:=$(shell R CMD config LDFLAGS)
CPICFLAGS:=$(shell R CMD config CPICFLAGS)

L += ${SOCKETLIB} -lm -lz

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
  [ -d jkweb/lib/${MATCHTYPE} ] || mkdir jkweb/lib/${MATCHTYPE}
	$(MAKE) -C jkweb/lib
	ln -s jkweb/lib/${MACHTYPE}/jkweb.a jkweb.a	

bigwig.so: prepare jkweb.a bigWig_R.c utils_R.c bwgExtra.h bwgExtra.c bw_query.h bw_query.c bw_base.h bw_base.c bigWig_R2.c
	@MAKEFLAGS=${MFLGSESC} R CMD SHLIB -o $@ bigWig_R.c utils_R.c bwgExtra.c  bw_base.c bw_query.c bigWig_R2.c jkweb.a ${L}

bigwiglib.a: jkweb.a bw_query.h bw_query.c bw_base.h bw_base.c bigwiglib.c
	$(CC) -fPIC -I./jkweb/inc ${COPT} ${CFLAGS} -c -o bigwiglib_base.o bw_base.c
	$(CC) -fPIC -I./jkweb/inc ${COPT} ${CFLAGS} -c -o bigwiglib_query.o bw_query.c
	$(CC) -fPIC -I./jkweb/inc ${COPT} ${CFLAGS} -c -o bigwiglib.o bigwiglib.c
	ar rcs $@ bigwiglib.o bigwiglib_base.o bigwiglib_query.o `find jkweb/ -name "*.o"`
	[ -d ../inst/libs ] || mkdir ../inst/libs
	cp bigwiglib.a ../inst/libs

bigwiglib.so: jkweb.a bw_query.h bw_query.c bw_base.h bw_base.c bigwiglib.c
	$(CC) -shared -fPIC -I./jkweb/inc ${COPT} ${CFLAGS} -o $@ bw_base.c bw_query.c bigwiglib.c jkweb.a ${L} ${LDFLAGS}

clean:
	rm -f bigwig.so bigwiglib.a bigwiglib.so *.o jkweb.a
	$(MAKE) -C jkweb/lib clean

tgtclean:
	rm -f bigwig.so bigwiglib.a bigwiglib.so *.o

.PHONY: prepare
