Calc Mappability
================

Scripts to identify positions in the target genome that cannot be mapped due to mapping ambiguities.

Requirements
============

1. Genome tools (URL: https://github.com/genometools/genometool).
2. Bedops package (URL: https://github.com/bedops/bedops).
3. Kent source utilities (URL: http://genome.ucsc.edu/FAQ/FAQlicense.html#license3).

Input
=====

1. Fasta and 2bit file representing the target genome.
2. The target read size. 

Output
======

1. A bed file representing the set of genomic coordinates where a read can not be unambiguously mapped.
2. A bigWig file, compatible with the bigWig package.
3. Several (large) temporary that can be removed after running.

Instructions
============

1. Download and install the necessary dependencies.
2. Modify the calcUnmap.sh and bed2bigWig.bsh scripts to update the target read size and genome location parameters.
3. Run the calcUnmap.sh script.
4. Run the bed2bigWig.bsh script.
