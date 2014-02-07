library(bigWig)

tmp = load.bigWig("test3.bigWig")

queryByStep_.bigWig <- function (bigWig, chrom, start, end, step, do.sum = FALSE, default.null = TRUE, defaultValue = 0) 
{
  op = "avg"
  if (do.sum)
    op = "sum"

  res = step.probeQuery.bigWig(bigWig, chrom, start, end, step, op = op, with.attributes = FALSE)

  if (all(is.na(res))) {
    if (default.null)
      return(NULL)

    return(rep(defaultValue, (end - start)%/%step))
  }

  # fill gaps with defaultValue
  res[is.na(res)] = defaultValue

  return(res)
}

vOld = queryByStep.bigWig(tmp, "chrU", 1, 101, 5)
vNew = queryByStep_.bigWig(tmp, "chrU", 1, 101, 5)
stopifnot(all.equal(vOld, vNew))

vOld = queryByStep.bigWig(tmp, "chrU", 1, 106, 15)
vNew = queryByStep_.bigWig(tmp, "chrU", 1, 106, 15)
stopifnot(all.equal(vOld, vNew))

vOld = queryByStep.bigWig(tmp, "chrU", 100, 104, 1)
vNew = queryByStep_.bigWig(tmp, "chrU", 100, 104, 1)
stopifnot(all.equal(vOld, vNew))

vOld = queryByStep.bigWig(tmp, "chrU", 1, 51, 5, do.sum=TRUE)
vNew = queryByStep_.bigWig(tmp, "chrU", 1, 51, 5, do.sum=TRUE)
stopifnot(all.equal(vOld, vNew))

vOld = queryByStep.bigWig(tmp, "chrU", 1, 61, 15)
vNew = queryByStep_.bigWig(tmp, "chrU", 1, 61, 15)
stopifnot(all.equal(vOld, vNew))

# example where bpQuery diverges in behaviour!
# step.bpQuery.bigWig(tmp, "chrU", 1, 61, 15, with.attributes=FALSE)
