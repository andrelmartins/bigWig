context("boundary tests")

test_that("region probe sum", {

  # start/end within probe
  # chrP 60 65 9
  res = region.probeQuery.bigWig(bwTest, "chrP", 61, 63, op = "sum")
  expect_equal(res, 9)

  # start/end covers probe
  res = region.probeQuery.bigWig(bwTest, "chrP", 58, 66, op = "sum")
  expect_equal(res, 9)

  # start @ edge of probe / end inside
  # chrP    10  20    10
  # chrP    20  30     5
  res = region.probeQuery.bigWig(bwTest, "chrP", 20, 25, op = "sum")
  expect_equal(res, 5)

  # start inside / end @ edge of probe
  res = region.probeQuery.bigWig(bwTest, "chrP", 15, 20, op = "sum")
  expect_equal(res, 10)

  # start inside / end over next probe
  res = region.probeQuery.bigWig(bwTest, "chrP", 15, 21, op = "sum")
  expect_equal(res, 15)

  # end after end of chrom
  # chrP 100
  # chrP    60  65     9
  res = region.probeQuery.bigWig(bwTest, "chrP", 60, 120, op = "sum")
  expect_equal(res, 9)

  # covers nothing (gap.value)
  res = region.probeQuery.bigWig(bwTest, "chrP", 70, 120, op = "sum")
  expect_true(is.na(res))

  # edges: start < 0
  expect_error(region.probeQuery.bigWig(bwTest, "chrP", -2, 10, op="sum"))

  # end > end of chrom
  # explicit gap value
  res = region.probeQuery.bigWig(bwTest, "chrP", 70, 120, op = "sum", gap.value = 2)
  expect_equal(res, 2)

  # end <= start
  expect_error(region.probeQuery.bigWig(bwTest, "chrP", 1, 1, op = "sum"))
  expect_error(region.probeQuery.bigWig(bwTest, "chrP", 3, 1, op = "sum"))
})


test_that("region base pair sum", {

  # start/end within probe
  # chrP 60 65 9
  res = region.bpQuery.bigWig(bwTest, "chrP", 61, 63, op = "sum")
  expect_equal(res, 9*(63 - 61))

  # start/end covers probe
  res = region.bpQuery.bigWig(bwTest, "chrP", 58, 66, op = "sum")
  expect_equal(res, 9*(65 - 60))

  # start @ edge of probe / end inside
  # chrP    10  20    10
  # chrP    20  30     5
  res = region.bpQuery.bigWig(bwTest, "chrP", 20, 25, op = "sum")
  expect_equal(res, 5*(25 - 20))

  # chrB 1 2 1
  # chrB 2 3 2
  res = region.bpQuery.bigWig(bwTest, "chrB", 2, 3, op = "sum")
  expect_equal(res, 2)

  # start inside / end @ edge of probe
  res = region.bpQuery.bigWig(bwTest, "chrP", 15, 20, op = "sum")
  expect_equal(res, 10 * (20 - 15))

  # start inside / end over next probe
  res = region.bpQuery.bigWig(bwTest, "chrP", 15, 21, op = "sum")
  expect_equal(res, 10*(20 - 15) + 5)

  # end after end of chrom
  # chrP 100
  # chrP    60  65     9
  res = region.bpQuery.bigWig(bwTest, "chrP", 60, 120, op = "sum")
  expect_equal(res, 9*(65 - 60))

  # covers nothing (gap.value)
  res = region.bpQuery.bigWig(bwTest, "chrP", 70, 120, op = "sum")
  expect_equal(res, 0) # default gap.value

  # edges: start < 0
  expect_error(region.bpQuery.bigWig(bwTest, "chrP", -2, 10, op="sum"))

  # end > end of chrom
  # explicit gap value
  res = region.bpQuery.bigWig(bwTest, "chrP", 70, 120, op = "sum", gap.value = 2)
  expect_equal(res, 2)

  # end <= start
  expect_error(region.bpQuery.bigWig(bwTest, "chrP", 1, 1, op = "sum"))
  expect_error(region.bpQuery.bigWig(bwTest, "chrP", 3, 1, op = "sum"))
})
