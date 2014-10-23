#
# compiler definitions for R
#

CC:=$(shell $(R_HOME)/bin/R CMD config CC)
MACHTYPE:=$(shell $(R_HOME)/bin$(R_ARCH_BIN)/R -e "cat('@', R.version[['arch']], '@', sep='')" | grep -o "^@.*@" | sed -e 's/^@\(.*\)@/\1/')
LDFLAGS:=$(shell $(R_HOME)/bin$(R_ARCH_BIN)/R CMD config LDFLAGS)
CPICFLAGS:=$(shell $(R_HOME)/bin$(R_ARCH_BIN)/R CMD config CPICFLAGS)
COPT:=
CFLAGS += $(shell $(R_HOME)/bin$(R_ARCH_BIN)/R CMD config CFLAGS) $(CPICFLAGS) -I${R_HOME}/include 
