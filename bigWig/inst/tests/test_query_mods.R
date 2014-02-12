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

# TODO: check for warning


# TODO: check 'end' attribute
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
# TODO
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

