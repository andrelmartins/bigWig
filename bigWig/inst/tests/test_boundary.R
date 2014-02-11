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


})
