suppressPackageStartupMessages({
  library(dplyr)
  library(testthat)
})

test_that("Overlap indices.", {
  
  # Case: 1:1 and 2:2 overlap.
  left <- data.frame(
    chr = c(1, 1, 1),
    start = c(1, 100, 1000),
    end = c(10, 110, 1010)
  )
  
  right <- data.frame(
    chr = c(1, 1, 1),
    start = c(5, 105, 1011),
    end = c(15, 108, 1015)
  )
  
  exp <- data.frame(
    left_idx = c(1, 2),
    right_idx = c(1, 2)
  )
  obs <- OverlapIndices(left, right)
  expect_equal(obs, exp)
  
  # Case: 1:2 and 1:3 overlap.
  left <- data.frame(
    chr = c(1, 2),
    start = c(1, 1),
    end = c(1000, 1000)
  )
  
  right <- data.frame(
    chr = c(0, 1, 1),
    start = c(1, 10, 100),
    end = c(1, 20, 200)
  )
  
  exp <- data.frame(
    left_idx = c(1, 1),
    right_idx = c(2, 3)
  )
  obs <- suppressWarnings(OverlapIndices(left, right))
  expect_equal(obs, exp)
  
  # Case: no overlaps.
  left <- data.frame(
    chr = c(1),
    start = c(1),
    end = c(10)
  )
  
  right <- data.frame(
    chr = c(1, 2),
    start = c(20, 1),
    end = c(30, 10)
  )
  
  exp <- data.frame(
    left_idx = integer(),
    right_idx = integer()
  )
  obs <- OverlapIndices(left, right)
  expect_equal(obs, exp)
  
})


# -----------------------------------------------------------------------------


test_that("Find overlaps.", {
  
  # Case: 1:1 and 2:2 overlap.
  left <- data.frame(
    chr = c(1, 1, 1),
    start = c(1, 100, 1000),
    end = c(10, 110, 1010)
  )
  
  right <- data.frame(
    chr = c(1, 1, 1),
    start = c(5, 105, 1011),
    end = c(15, 108, 1015)
  )
  
  exp <- data.frame(
    left_chr = c(1, 1, 1),
    left_start = c(1, 100, 1000),
    left_end = c(10, 110, 1010),
    left_idx = c(1, 2, 3),
    right_idx = c(1, 2, NA),
    right_chr = c(1, 1, NA),
    right_start = c(5, 105, NA),
    right_end = c(15, 108, NA),
    any_overlaps = c(TRUE, TRUE, FALSE)
  )
  obs <- FindOverlaps(left, right)
  expect_equal(obs, exp)
  
  # Reversing arguments to FindOverlaps.
  exp <- data.frame(
    left_chr = c(1, 1, 1),
    left_start = c(5, 105, 1011),
    left_end = c(15, 108, 1015),
    left_idx = c(1, 2, 3),
    right_idx = c(1, 2, NA),
    right_chr = c(1, 1, NA),
    right_start = c(1, 100, NA),
    right_end = c(10, 110, NA),
    any_overlaps = c(TRUE, TRUE, FALSE)
  )
  obs <- FindOverlaps(right, left)
  expect_equal(obs, exp)
  
  # Case: no overlaps.
  left <- data.frame(
    chr = c(1),
    start = c(1),
    end = c(10)
  )
  
  right <- data.frame(
    chr = c(1, 2),
    start = c(20, 1),
    end = c(30, 10)
  )
  
  exp <- data.frame(
    left_chr = 1,
    left_start = 1,
    left_end = 10,
    left_idx = 1,
    right_idx = as.integer(NA),
    right_chr = as.integer(NA),
    right_start = as.integer(NA),
    right_end = as.integer(NA),
    any_overlaps = FALSE
  )
  obs <- FindOverlaps(left, right)
  expect_equal(obs, exp)
  
})


# -----------------------------------------------------------------------------


test_that("Calculate overlap statistics.", {
  
  # Case: 1:1, 1:2, 3:3.
  left <- data.frame(
    chr = c(1, 1, 1),
    start = c(1, 200, 400),
    end = c(100, 300, 500)
  )
  
  right <- data.frame(
    chr = c(1, 1, 1),
    start = c(20, 50, 440),
    end = c(25, 55, 450)
  )
  
  out <- CalcOverlapStats(left, right)
  expect_equal(out@left_overlapping, c(1, 3))
  expect_equal(out@right_overlapped, c(1, 2, 3))
  
})
