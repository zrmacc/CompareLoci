#' Overlap Indices
#'
#' Identifies which loci from \code{left} overlap with which
#' loci from \code{right}.
#'
#' @param left Data.frame containing the first set of loci.
#' @param right Data.frame containing the second set of loci.
#' @param left_chr_name Name of column containing chromosome.
#' @param left_start_name Name of column containing locus start position.
#' @param left_end_name Name of column containing locus end position.
#' @param right_chr_name Name of column containing chromosome.
#' @param right_start_name Name of column containing locus start position.
#' @param right_end_name Name of column containing locus end position.
#' @return Data.frame.
#' @export 

OverlapIndices <- function(
    left,
    right,
    left_chr_name = "chr",
    left_start_name = "start",
    left_end_name = "end",
    right_chr_name = "chr",
    right_start_name = "start",
    right_end_name = "end"
) {
  gr1 <- DF2GR(left, left_chr_name, left_start_name, left_end_name)
  gr2 <- DF2GR(right, right_chr_name, right_start_name, right_end_name)
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
#' @param left_chr_name Name of column containing chromosome.
#' @param left_start_name Name of column containing locus start position.
#' @param left_end_name Name of column containing locus end position.
#' @param right_chr_name Name of column containing chromosome.
#' @param right_start_name Name of column containing locus start position.
#' @param right_end_name Name of column containing locus end position.
#' @return Data.frame. 
#' @importFrom dplyr all_of "%>%"
#' @export 

FindOverlaps <- function(
    left,
    right,
    left_chr_name = "chr",
    left_start_name = "start",
    left_end_name = "end",
    right_chr_name = "chr",
    right_start_name = "start",
    right_end_name = "end"
) {
  if (!methods::is(left, "data.frame") & !methods::is(right, "data.frame")) {
    stop("Inputs should be data frames.")
  }
  
  # Rename chr, start, end.
  left <- left %>%
    dplyr::rename(
      chr = all_of(left_chr_name),
      start = all_of(left_start_name),
      end = all_of(left_end_name)
    )
  right <- right %>%
    dplyr::rename(
      chr = all_of(right_chr_name),
      start = all_of(right_start_name),
      end = all_of(right_end_name)
    )
  
  # Overlaps.
  overlap_indices <- OverlapIndices(left, right)
  
  # Prefix right.
  right$idx <- seq_len(nrow(right))
  colnames(right) <- paste0("right_", colnames(right))
  
  # Inner join overlaps with right.
  overlap_indices <- overlap_indices %>%
    dplyr::inner_join(right, by = "right_idx")
  
  # Prefix left.
  left$idx <- seq_len(nrow(left))
  colnames(left) <- paste0("left_", colnames(left))
  
  # Left join left loci with overlaps.
  right_idx <- NULL
  out <- left %>%
    dplyr::left_join(overlap_indices, by = "left_idx") %>%
    dplyr::mutate(
      any_overlaps = !is.na(right_idx)
    ) %>%
    dplyr::relocate(
      "left_idx", "left_chr", "left_start", "left_end",
      "right_idx", "right_chr", "right_start", "right_end"
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
#' @param left_chr_name Name of column containing chromosome.
#' @param left_start_name Name of column containing locus start position.
#' @param left_end_name Name of column containing locus end position.
#' @param right_chr_name Name of column containing chromosome.
#' @param right_start_name Name of column containing locus start position.
#' @param right_end_name Name of column containing locus end position.
#' @return Data.frame.
#' @importFrom dplyr all_of "%>%"
#' @export 

CalcOverlapStats <- function(
    left,
    right,
    left_chr_name = "chr",
    left_start_name = "start",
    left_end_name = "end",
    right_chr_name = "chr",
    right_start_name = "start",
    right_end_name = "end"
) {
  
  # Rename chr, start, end.
  left <- left %>%
    dplyr::rename(
      chr = all_of(left_chr_name),
      start = all_of(left_start_name),
      end = all_of(left_end_name)
    )
  right <- right %>%
    dplyr::rename(
      chr = all_of(right_chr_name),
      start = all_of(right_start_name),
      end = all_of(right_end_name)
    )
  
  # Find overlaps.
  overlaps <- FindOverlaps(left, right)
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