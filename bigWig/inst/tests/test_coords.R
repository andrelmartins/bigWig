context("coordinate checks")

test_that("chrom missing", {

  bed = data.frame("chrBogus", 10, 20)
  bed6 = data.frame("chrBogus", 10, 20, "N", 0, "+")

  expect_error(region.bpQuery.bigWig(bwTest, "chrBogus", 10, 20))
  expect_error(bed.region.bpQuery.bigWig(bwTest, bed))
  expect_error(bed6.region.bpQuery.bigWig(bwTest, bwTest, bed6))
  expect_error(step.bpQuery.bigWig(bwTest, "chrBogus", 10, 20, 1))
  expect_error(step.bpQuery.bigWig(bwTest, "chrBogus", NULL, NULL, 1))
  expect_error(bed.step.bpQuery.bigWig(bwTest, bed, 1))
  expect_error(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 1))

  expect_error(region.probeQuery.bigWig(bwTest, "chrBogus", 10, 20))
  expect_error(bed.region.probeQuery.bigWig(bwTest, bed))
  expect_error(bed6.region.probeQuery.bigWig(bwTest, bwTest, bed6))
  expect_error(step.probeQuery.bigWig(bwTest, "chrBogus", 10, 20, 1))
  expect_error(step.probeQuery.bigWig(bwTest, "chrBogus", NULL, NULL, 1))
  expect_error(bed.step.probeQuery.bigWig(bwTest, bed, 1))
  expect_error(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 1))

  #
  # same but now with alternative bw spec
  basePath = system.file("extdata", package="bigWig")
  prefix = file.path(basePath, "test3_")
  suffix = ".bigWig"
  bw = c(prefix, suffix)

  expect_error(region.bpQuery.bigWig(bw, "chrBogus", 10, 20))
  expect_error(bed.region.bpQuery.bigWig(bw, bed))
  expect_error(bed6.region.bpQuery.bigWig(bw, bw, bed6))
  expect_error(step.bpQuery.bigWig(bw, "chrBogus", 10, 20, 1))
  expect_error(step.bpQuery.bigWig(bw, "chrBogus", NULL, NULL, 1))
  expect_error(bed.step.bpQuery.bigWig(bw, bed, 1))
  expect_error(bed6.step.bpQuery.bigWig(bw, bw, bed6, 1))

  expect_error(region.probeQuery.bigWig(bw, "chrBogus", 10, 20))
  expect_error(bed.region.probeQuery.bigWig(bw, bed))
  expect_error(bed6.region.probeQuery.bigWig(bw, bw, bed6))
  expect_error(step.probeQuery.bigWig(bw, "chrBogus", 10, 20, 1))
  expect_error(step.probeQuery.bigWig(bw, "chrBogus", NULL, NULL, 1))
  expect_error(bed.step.probeQuery.bigWig(bw, bed, 1))
  expect_error(bed6.step.probeQuery.bigWig(bw, bw, bed6, 1))
})

test_that("end > start", {

  bed = data.frame("chrB", c(10, 11), c(9, 11))
  bed6 = data.frame("chrB", c(10, 11), c(9, 11), "N", 0, "+")

  expect_error(region.bpQuery.bigWig(bwTest, "chrB", 10, 9))
  expect_error(region.bpQuery.bigWig(bwTest, "chrB", 11, 11))
  expect_error(bed.region.bpQuery.bigWig(bwTest, bed))
  expect_error(bed6.region.bpQuery.bigWig(bwTest, bwTest, bed6))
  expect_error(step.bpQuery.bigWig(bwTest, "chrB", 10, 9, 1))
  expect_error(step.bpQuery.bigWig(bwTest, "chrB", 11, 11, 1))
  expect_error(bed.step.bpQuery.bigWig(bwTest, bed, 1))
  expect_error(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 1))

  expect_error(region.probeQuery.bigWig(bwTest, "chrB", 10, 9))
  expect_error(region.probeQuery.bigWig(bwTest, "chrB", 11, 11))
  expect_error(bed.region.probeQuery.bigWig(bwTest, bed))
  expect_error(bed6.region.probeQuery.bigWig(bwTest, bwTest, bed6))
  expect_error(step.probeQuery.bigWig(bwTest, "chrB", 10, 9, 1))
  expect_error(step.probeQuery.bigWig(bwTest, "chrB", 11, 11, 1))
  expect_error(bed.step.probeQuery.bigWig(bwTest, bed, 1))
  expect_error(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 1))
})

test_that("start >= 0", {

  bed = data.frame("chrB", -2, 10)
  bed6 = data.frame("chrB", -1, 10, "N", 0, "+")

  expect_error(region.bpQuery.bigWig(bwTest, "chrB", -2, 9))
  expect_error(bed.region.bpQuery.bigWig(bwTest, bed))
  expect_error(bed6.region.bpQuery.bigWig(bwTest, bwTest, bed6))
  expect_error(step.bpQuery.bigWig(bwTest, "chrB", -2, 9, 1))
  expect_error(bed.step.bpQuery.bigWig(bwTest, bed, 1))
  expect_error(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 1))

  expect_error(region.probeQuery.bigWig(bwTest, "chrB", -2, 9))
  expect_error(bed.region.probeQuery.bigWig(bwTest, bed))
  expect_error(bed6.region.probeQuery.bigWig(bwTest, bwTest, bed6))
  expect_error(step.probeQuery.bigWig(bwTest, "chrB", 10, -2, 1))
  expect_error(bed.step.probeQuery.bigWig(bwTest, bed, 1))
  expect_error(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 1))

})

test_that("invalid strand", {

  bed = data.frame("chrB", 10, 20)
  bed6 = data.frame("chrB", 10, 20, "N", 0, "?")

  expect_error(region.bpQuery.bigWig(bwTest, "chrB", 10, 20, strand = "?"))
  expect_error(bed.region.bpQuery.bigWig(bwTest, bed, strand = "?"))
  expect_error(bed6.region.bpQuery.bigWig(bwTest, bwTest, bed6))
  expect_error(step.bpQuery.bigWig(bwTest, "chrB", 10, 20, 1, strand = "?"))
  expect_error(step.bpQuery.bigWig(bwTest, "chrB", NULL, NULL, 1, strand = "?"))
  expect_error(bed.step.bpQuery.bigWig(bwTest, bed, 1, strand = "?"))
  expect_error(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 1))

  expect_error(bed6.region.probeQuery.bigWig(bwTest, bwTest, bed6))
  expect_error(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 1))
})

test_that("no strand in BED6 queries", {
  bed5 = data.frame("chrB", 10, 20, "N", 0)

  expect_error(bed6.region.bpQuery.bigWig(bwTest, bwTest, bed5))
  expect_error(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed5, 1))

  expect_error(bed6.region.probeQuery.bigWig(bwTest, bwTest, bed5))
  expect_error(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed5, 1))
})
