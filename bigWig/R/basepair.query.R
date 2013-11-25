#
# Query functions for "base pair" mode
#

valid.strand <- function(strand) {
  !is.na(strand) && (strand == "+" || strand == "-")
}

region.bpQuery.bigWig <- function(bw, chrom, start, end, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL) {
  if (!is.null(bwMap) && !valid.strand(strand))
    error("strand is required when using mappability information")

}

bed.region.bpQuery.bigWig <- function(bw, bed, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL) {
  if (!is.null(bwMap) && !valid.strand(strand))
    error("strand is required when using mappability information")
  
}

bed6.region.bpQuery.bigWig <- function(bw, bed, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL) {
  
}

# note: start, end are optional here (use NULL for both to get the entire choromosome)
step.bpQuery.bigWig <- function(bw, chrom, start, end, step, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL, with.attributes = TRUE) {
  if (!is.null(bwMap) && !valid.strand(strand))
    error("strand is required when using mappability information")
  
}

bed.step.bpQuery.bigWig <- function(bw, bed, step, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL, with.attributes = FALSE, as.matrix = FALSE) {
  if (!is.null(bwMap) && !valid.strand(strand))
    error("strand is required when using mappability information")
  
}

bed6.step.bpQuery.bigWig <- function(bw.plus, bw.minus, bed6, step, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL, with.attributes = FALSE, as.matrix = FALSE) {
  
}
