#!/usr/bin/sh
#
# This script calculates the set of unique positions that are identical in a target genome.  
# 
# Requirements (In the current path): 
# * The genometools package (tested with 1.5.1).
# * The bedops package.
# 
# Input: 
# * A mappability bed file created by calcUnmap.sh.
# * A 2bit formatted genome file.
#
# Output:
# * A bigWig of the positions in the target genome that are repetitive
#   at the target read size. 

SCRATCH="."
GENOME="genomes/hg19/hg19.fa.gz"
mapsize=30

## Create a representation of the target genome using suffix trees.
gt suffixerator -dna -pl -tis -suf -lcp -v -parts 4 -db $GENOME -indexname $SCRATCH/reads

## Make an index of the genome at a given read size.
gt tallymer mkindex -mersize $mapsize -minocc 2  -indexname $SCRATCH/$mapsize\mers  -counts -pl -esa $SCRATCH/reads

## Write the positions for sequences that occur two or more times into a useful text based file.
gt tallymer search -output qseqnum qpos counts sequence -strand fp -tyr $SCRATCH/$mapsize\mers -q $GENOME | gzip > $SCRATCH/$mapsize\mers.gtTxt.gz

## Transfer to a (MUCH) more compact bed format.
SEQNAMES=`zcat $GENOME | grep ">" | sed "s/^>//g"`
zcat $SCRATCH/$mapsize\mers.gtTxt.gz | perl tallymer2bed.pl $SEQNAMES | bedops -m  - | gzip > $mapsize\mers.unmap.bed.gz


