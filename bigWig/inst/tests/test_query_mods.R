context("query modifiers")

test_that("step >= 1", {

bed = data.frame("chrB", 10, 20)
bed6 = data.frame("chrB", 10, 20, "N", 0, "+")

for (step in c(-1, 0)) {
  expect_error(step.bpQuery.bigWig(bwTest, "chrB", 10, 20, step))
  expect_error(step.bpQuery.bigWig(bwTest, "chrB", NULL, NULL, step))
  expect_error(bed.step.bpQuery.bigWig(bwTest, bed, step))
  expect_error(bed6.step.bpQuery.bigWig(bwTest, bed6, step))

  expect_error(step.probeQuery.bigWig(bwTest, "chrB", 10, 20, step))
  expect_error(step.probeQuery.bigWig(bwTest, "chrB", NULL, NULL, step))
  expect_error(bed.step.probeQuery.bigWig(bwTest, bed, step))
  expect_error(bed6.step.probeQuery.bigWig(bwTest, bed6, step))
}
})

test_that("query not multiple of step size", {

# TODO: check for warning


# TODO: check 'end' attribute
})


test_that("abs.value", {
# TODO
})

test_that("strand reverse result", {
# TODO
})

test_that("as.matrix", {
# TODO: result is matrix or not

# TODO: result attributes
})

