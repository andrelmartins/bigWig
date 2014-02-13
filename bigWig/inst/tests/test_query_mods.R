context("query modifiers")

test_that("step >= 1", {

  bed = data.frame("chrB", 10, 20)
  bed6 = data.frame("chrB", 10, 20, "N", 0, "+")

  for (step in c(-1, 0)) {
    expect_error(step.bpQuery.bigWig(bwTest, "chrB", 10, 20, step))
    expect_error(step.bpQuery.bigWig(bwTest, "chrB", NULL, NULL, step))
    expect_error(bed.step.bpQuery.bigWig(bwTest, bed, step))
    expect_error(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, step))

    expect_error(step.probeQuery.bigWig(bwTest, "chrB", 10, 20, step))
    expect_error(step.probeQuery.bigWig(bwTest, "chrB", NULL, NULL, step))
    expect_error(bed.step.probeQuery.bigWig(bwTest, bed, step))
    expect_error(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, step))
  }
})

test_that("query not multiple of step size", {

  bed = data.frame("chrB", 10, 20)
  bed6 = data.frame("chrB", 10, 20, "N", 0, "+")
  
  # check for warning & 'end' attribute correctness
  expect_warning((v1 = step.bpQuery.bigWig(bwTest, "chrB", 10, 20, 3)), "query region is not an integer multiple of step")
  expect_equal(attr(v1, "end"), 10 + length(v1)*3)
  
  expect_warning((v1 = bed.step.bpQuery.bigWig(bwTest, bed, 3, with.attributes = TRUE)), "query region is not an integer multiple of step")
  for(res in v1)
    expect_equal(attr(res, "end"), 10 + length(res)*3)
  
  expect_warning((v1 = bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 3, with.attributes = TRUE)), "query region is not an integer multiple of step")
  for(res in v1)
    expect_equal(attr(res, "end"), 10 + length(res)*3)
  
  #
  expect_warning((v1 = step.probeQuery.bigWig(bwTest, "chrB", 10, 20, 3)), "query region is not an integer multiple of step")
  expect_equal(attr(v1, "end"), 10 + length(v1)*3)
  
  expect_warning((v1 = bed.step.probeQuery.bigWig(bwTest, bed, 3, with.attributes = TRUE)), "query region is not an integer multiple of step")
  for(res in v1)
    expect_equal(attr(res, "end"), 10 + length(res)*3)
  
  expect_warning((v1 = bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 3, with.attributes = TRUE)), "query region is not an integer multiple of step")
  for(res in v1)
    expect_equal(attr(res, "end"), 10 + length(res)*3)
})


test_that("abs.value", {

  bed = data.frame(c("chrB", "chrBm"), 10, 20)
  bed6 = data.frame(c("chrB", "chrBm"), 10, 20, "N", 0, "+")

  expect_true(region.bpQuery.bigWig(bwTest, "chrBm", 10, 20) < 0)
  expect_true(region.bpQuery.bigWig(bwTest, "chrBm", 10, 20, abs.value = TRUE) >= 0)

  expect_equal(bed.region.bpQuery.bigWig(bwTest, bed) < 0,  c(FALSE, TRUE))
  expect_equal(bed.region.bpQuery.bigWig(bwTest, bed, abs.value = TRUE) < 0, c(FALSE, FALSE))

  expect_equal(bed6.region.bpQuery.bigWig(bwTest, bwTest, bed6) < 0, c(FALSE, TRUE))
  expect_equal(bed6.region.bpQuery.bigWig(bwTest, bwTest, bed6, abs.value = TRUE) < 0, c(FALSE, FALSE))

  expect_true(all(step.bpQuery.bigWig(bwTest, "chrBm", 10, 20, 1) <= 0))
  expect_true(all(step.bpQuery.bigWig(bwTest, "chrBm", 10, 20, 1, abs.value = TRUE) >= 0))

  expect_true(all(bed.step.bpQuery.bigWig(bwTest, bed, 1, as.matrix = TRUE)[2,] <= 0))
  expect_true(all(bed.step.bpQuery.bigWig(bwTest, bed, 1, as.matrix = TRUE, abs.value = TRUE) >= 0))

  expect_true(all(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 1, as.matrix = TRUE)[2,] <= 0))
  expect_true(all(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 1, as.matrix = TRUE, abs.value = TRUE) >= 0))

  #
  # now with probe functions
  expect_true(region.probeQuery.bigWig(bwTest, "chrBm", 10, 20) < 0)
  expect_true(region.probeQuery.bigWig(bwTest, "chrBm", 10, 20, abs.value = TRUE) >= 0)

  expect_equal(bed.region.probeQuery.bigWig(bwTest, bed) < 0,  c(FALSE, TRUE))
  expect_equal(bed.region.probeQuery.bigWig(bwTest, bed, abs.value = TRUE) < 0, c(FALSE, FALSE))

  expect_equal(bed6.region.probeQuery.bigWig(bwTest, bwTest, bed6) < 0, c(FALSE, TRUE))
  expect_equal(bed6.region.probeQuery.bigWig(bwTest, bwTest, bed6, abs.value = TRUE) < 0, c(FALSE, FALSE))

  expect_true(all(step.probeQuery.bigWig(bwTest, "chrBm", 10, 20, 1) < 0, na.rm = TRUE))
  expect_true(all(step.probeQuery.bigWig(bwTest, "chrBm", 10, 20, 1, abs.value = TRUE) > 0, na.rm = TRUE))

  expect_true(all(bed.step.probeQuery.bigWig(bwTest, bed, 1, as.matrix = TRUE)[2,] < 0, na.rm = TRUE))
  expect_true(all(bed.step.probeQuery.bigWig(bwTest, bed, 1, as.matrix = TRUE, abs.value = TRUE) > 0, na.rm = TRUE))

  expect_true(all(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 1, as.matrix = TRUE)[2,] < 0, na.rm = TRUE))
  expect_true(all(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 1, as.matrix = TRUE, abs.value = TRUE) > 0, na.rm = TRUE))
})

test_that("strand reverse result", {
  bed = data.frame(c("chrB", "chrBm"), 10, 20)
  bed6 = data.frame(c("chrB", "chrBm"), 10, 20, "N", 0, "+")
  bed6.m = data.frame(c("chrB", "chrBm"), 10, 20, "N", 0, "-")

  # in step queries, - strand does *not* reverse the result
  v1 = step.bpQuery.bigWig(bwTest, "chrB", 10, 20, 1)
  v2 = step.bpQuery.bigWig(bwTest, "chrB", 10, 20, 1, strand = '-')
  expect_equal(v1, v2)
  
  # full chrom query
  v1 = step.bpQuery.bigWig(bwTest, "chrB", NULL, NULL, 1)
  v2 = step.bpQuery.bigWig(bwTest, "chrB", NULL, NULL, 1, strand = '-')
  expect_equal(v1, v2)
  
  v1 = bed.step.bpQuery.bigWig(bwTest, bed, 1)
  v2 = bed.step.bpQuery.bigWig(bwTest, bed, 1, strand = '-')
  for (i in 1:(dim(bed)[1]))
    expect_equal(v1[[i]], v2[[i]])
  
  # but bed6 do with follow.strand = TRUE
  v1 = bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 1)
  v2 = bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6.m, 1)
  v3 = bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6.m, 1, follow.strand = TRUE)
  for (i in 1:(dim(bed)[1])) {
    expect_equal(v1[[i]], v2[[i]])
    expect_equal(as.vector(v1[[i]]), rev(v3[[i]]))
  }
  
  # now for probes (only bed6 has strand, same as above)
  v1 = bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 1)
  v2 = bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6.m, 1)
  v3 = bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6.m, 1, follow.strand = TRUE)
  for (i in 1:(dim(bed)[1])) {
    expect_equal(v1[[i]], v2[[i]])
    expect_equal(as.vector(v1[[i]]), rev(v3[[i]]))
  }
})

test_that("as.matrix", {
  #
  # this only applies to bed.step.* and bed6.step.* operations

  bed = data.frame(c("chrB", "chrBm"), 10, 20)
  bed6 = data.frame(c("chrB", "chrBm"), 10, 20, "N", 0, "+")

  expect_is(bed.step.bpQuery.bigWig(bwTest, bed, 1), "list")
  expect_is(bed.step.bpQuery.bigWig(bwTest, bed, 1, as.matrix = TRUE), "matrix")

  expect_is(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 1), "list")
  expect_is(bed6.step.bpQuery.bigWig(bwTest, bwTest, bed6, 1, as.matrix = TRUE), "matrix")

  expect_is(bed.step.probeQuery.bigWig(bwTest, bed, 1), "list")
  expect_is(bed.step.probeQuery.bigWig(bwTest, bed, 1, as.matrix = TRUE), "matrix")

  expect_is(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 1), "list")
  expect_is(bed6.step.probeQuery.bigWig(bwTest, bwTest, bed6, 1, as.matrix = TRUE), "matrix")
})

