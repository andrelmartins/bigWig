library(bigWig)

fpath = system.file("extdata", "test3.bigWig", package="bigWig")
tmp = load.bigWig(fpath)

chromStepSum_.bigWig <- function (bigwig, chrom, step, defaultValue) {
  step.bpQuery.bigWig(bigwig, chrom, NULL, NULL, step, gap.value = defaultValue, with.attributes=FALSE)
}

vOld = chromStepSum.bigWig(tmp, 'chrU', 5, 0)
vNew = chromStepSum_.bigWig(tmp, 'chrU', 5, 0)
stopifnot(all.equal(vOld, vNew))

vOld = chromStepSum.bigWig(tmp, 'chrU', 1, 1)
vNew = chromStepSum_.bigWig(tmp, 'chrU', 1, 1)
stopifnot(all.equal(vOld, vNew))
