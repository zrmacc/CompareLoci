suppressPackageStartupMessages({
  library(dplyr)
  library(testthat)
})

test_that("Data frame to genomic ranges.", {
  
  data <- data.frame(
    chr = c(1, 2, 3),
    start = c(101, 201, 301),
    end = c(110, 220, 330),
    covar = c(1.1, 2.2, 3.3)
  )
  
  gr <- DF2GR(data)
  expect_true(is(gr, "GRanges"))
  
  data <- data %>%
    dplyr::rename(
      CHR = chr,
      START = start,
      END = end
    )
  
  # With renaming.
  gr <- DF2GR(
    data,
    chr_name = "CHR",
    start_name = "START",
    end_name = "END"
  )
  expect_true(is(gr, "GRanges")) 
  
  # Incorrect input type.
  expect_error(DF2GR(gr))
  
})


# -----------------------------------------------------------------------------


test_that("Genomic ranges to data frame.", {
  
  data <- data.frame(
    chr = c(1, 2, 3),
    start = c(101, 201, 301),
    end = c(110, 220, 330),
    covar = c(1.1, 2.2, 3.3)
  )
  
  gr <- DF2GR(data)
  df <- GR2DF(gr)
  expect_true(is(df, "data.frame")) 
  
  df <- GR2DF(
    gr,
    chr_name = "CHR",
    start_name = "START",
    end_name = "END"
  )
  
  exp_colnames <- c("CHR", "START", "END", "covar")
  obs_colnames <- colnames(df)
  expect_true(all(exp_colnames == obs_colnames))
  
  expect_error(GR2DF(df))
  
})