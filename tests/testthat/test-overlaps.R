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
    left_idx = c(1, 2, 3),
    right_idx = c(1, 2, NA),
    any_overlaps = c(TRUE, TRUE, FALSE)
  )
  exp <- cbind(left, exp)
  obs <- FindOverlaps(left, right)
  expect_equal(obs, exp)
  
  # Non-standard column names.
  colnames(left) <- toupper(colnames(left))
  colnames(right) <- toupper(colnames(right))
  obs <- FindOverlaps(
    left,
    right,
    chr_name = "CHR",
    start_name = "START",
    end_name = "END"
  )
  colnames(exp)[1:3] <- toupper(colnames(exp)[1:3])
  expect_equal(obs, exp)
  
  # Reversing arguments to FindOverlaps.
  colnames(left) <- tolower(colnames(left))
  colnames(right) <- tolower(colnames(right))
  exp <- data.frame(
    left_idx = c(1, 2, 3),
    right_idx = c(1, 2, NA),
    any_overlaps = c(TRUE, TRUE, FALSE)
  )
  exp <- cbind(right, exp)
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
    left_idx = 1,
    right_idx = as.integer(NA),
    any_overlaps = FALSE
  )
  exp <- cbind(left, exp)
  obs <- FindOverlaps(left, right)
  expect_equal(obs, exp)
  
  # Case: 1:1 and 1:2. 
  left <- data.frame(
    chr = c(1, 1),
    start = c(1, 20),
    end = c(10, 30)
  )
  
  right <- data.frame(
    chr = c(1, 1),
    start = c(1, 3),
    end = c(2, 4)
  )
  
  exp <- data.frame(
    chr = c(1, 1, 1),
    start = c(1, 1, 20),
    end = c(10, 10, 30),
    left_idx = c(1, 1, 2),
    right_idx = c(1, 2, NA),
    any_overlaps = c(TRUE, TRUE, FALSE)
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