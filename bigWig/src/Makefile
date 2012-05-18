ifeq ($(KENTHOME),)
	KENTHOME :=/home/browser/kent-trunk
endif

include $(KENTHOME)/src/inc/common.mk

ifeq (${MACHTYPE},)
	MACHTYPE=x86_64
endif

L += ${SOCKETLIB} -lm -lz
MYLIBDIR = ./lib/${MACHTYPE}
MYLIBS = ${MYLIBDIR}/jkweb.a

# build our MAKEFLAGS
# R needs spaces to be escaped ...
#
MFLGS := "CFLAGS=-I./inc ${COPT} ${CFLAGS}"
# define a space
empty:=
space:= $(empty) $(empty)
MFLGSESC := $(subst $(space),\$(space),$(MFLGS))

all: bigwig.so

bigwig.so: bigWig_R.c
	ln -s $(KENTHOME)/src/inc
	ln -s $(KENTHOME)/src/lib
	MAKEFLAGS=${MFLGSESC} R CMD SHLIB -o $@ bigWig_R.c ${MYLIBS} ${L}
	rm inc lib

clean:
	rm -f bigwig.so *.o lib inc
