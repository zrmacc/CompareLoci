#' Overlap Indices
#'
#' Identifies which loci from \code{left} overlap with which
#' loci from \code{right}.
#'
#' @param left Data.frame containing the first set of loci.
#' @param right Data.frame containing the second set of loci.
#' @param chr_name Name of column containing chromosome.
#' @param start_name Name of column containing locus start position.
#' @param end_name Name of column containing locus end position.
#' @return Data.frame.
#' @export 

OverlapIndices <- function(
    left,
    right,
    chr_name = "chr",
    start_name = "start",
    end_name = "end"
) {
  gr1 <- DF2GR(left, chr_name, start_name, end_name)
  gr2 <- DF2GR(right, chr_name, start_name, end_name)
  indices <- as.data.frame(GenomicRanges::findOverlaps(gr1, gr2))
  colnames(indices) <- c("left_idx", "right_idx")
  return(indices)
}


# -----------------------------------------------------------------------------


#' Find Overlaps
#'
#' Identifies which loci in \code{left} overlap with which loci in \code{right}.
#' Note this function is asymmetric; \code{FindOverlaps(left, right)} will
#' generally differ from \code{FindOverlaps(right, left)}.
#'
#' @param left Data.frame containing the first set of loci.
#' @param right Data.frame containing the second set of loci.
#' @param chr_name Name of column containing chromosome.
#' @param start_name Name of column containing locus start position.
#' @param end_name Name of column containing locus end position.
#' @return Data.frame. 
#' @importFrom dplyr "%>%"
#' @export 

FindOverlaps <- function(
    left,
    right,
    chr_name = "chr",
    start_name = "start",
    end_name = "end"
) {
  if (!methods::is(left, "data.frame") & !methods::is(right, "data.frame")) {
    stop("Inputs should be data frames.")
  }
  out <- left
  out$left_idx <- seq_len(nrow(out))
  overlap_indices <- OverlapIndices(
    left, right, chr_name, start_name, end_name)
  right_idx <- NULL
  out <- out %>%
    dplyr::left_join(overlap_indices, by = "left_idx") %>%
    dplyr::mutate(
      any_overlaps = !is.na(right_idx)
    )
  return(out)
}


# -----------------------------------------------------------------------------


#' Calculate Overlap Statistics
#'
#' Calculates overlap statistics comparing loci on the \code{left} with loci on
#' the \code{right}. Note this function is asymmetric;
#' \code{CalcOverlapStats(left, right)} will generally differ from
#' \code{CalcOverlapStats(right, left)}.
#' 
#' @param left Data.frame containing the first set of loci.
#' @param right Data.frame containing the second set of loci.
#' @param chr_name Name of column containing chromosome.
#' @param start_name Name of column containing locus start position.
#' @param end_name Name of column containing locus end position.
#' @return Data.frame.
#' @export 

CalcOverlapStats <- function(
    left,
    right,
    chr_name = "chr",
    start_name = "start",
    end_name = "end"
) {
  
  overlaps <- FindOverlaps(left, right, chr_name, start_name, end_name)
  n_left <- nrow(left)
  n_right <- nrow(right)
  
  # Number of loci on the right overlapped by a locus on the left.
  right_overlapped <- setdiff(unique(overlaps$right_idx), NA)
  
  # Number of loci on the left that overlap at least 1 locus on the right.
  left_overlapping <- unique(overlaps$left_idx[overlaps$any_overlaps])
  
  out <- methods::new(
    Class = "overlap",
    left = left,
    right = right,
    overlaps = overlaps,
    left_overlapping = left_overlapping,
    right_overlapped = right_overlapped,
    n_left = n_left,
    n_right = n_right
  )
  
}