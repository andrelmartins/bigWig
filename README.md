bigWig
======

R interface to query UCSC bigWig files

Compilation
===========

1. Download the UCSC browser source (http://genome.ucsc.edu/FAQ/FAQlicense.html#license3)
2. Test that it compiles properly

  Look at the README file in the UCSC browser sources and follow up to step (6)

3. Edit bigWig/src/Makefile and change the KENTHOME to point to where you placed the UCSC browser source
4. Compile the package

        R CMD build --binary bigWig

    or

        R CMD install bigWig
