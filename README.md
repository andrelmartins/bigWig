bigWig
======
[![Build Status](https://travis-ci.org/andrelmartins/bigWig.svg?branch=master)](https://travis-ci.org/andrelmartins/bigWig)

R interface to query UCSC bigWig files

Compilation
===========

Simply do:

				R CMD INSTALL bigWig

By default, the CPU architecture used is the one used to build R (R.version$arch), unless the environment variable MACHTYPE is defined. If you need to specify a different architecture, export the environment variable MACHTYPE before compiling.
The default MACHTYPE is often a long string: "i386-redhat-linux-gnu" which will not function correctly in this build environment.
It needs to be something simple such as one of:

        i386 i686 sparc alpha x86_64 ppc etc ...

with no other alpha characters such as: -
To determine what your system reports itself as, try the uname options:  'uname -m' or 'uname -p' or 'uname -a' on your command line.  If necessary set this environment variable.

Do this under the bash shell as so:

       MACHTYPE=something
       export MACHTYPE

or under tcsh as so:

       setenv MACHTYPE something

Copyrights
==========

The code in src/jkweb is Copyright (c) 2000-2002 Jim Kent (c) 2003-2014 Regents of the University of California.  All other code is Copyright (c) 2012-2014 Cornell University.

See bigWig/DESCRIPTION, bigWig/LICENSE, bigWig/src/jkweb/README, bigWig/src/jkweb/LICENSE for license details.
