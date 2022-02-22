#' Convert from Data Frame to Genomic Ranges
#'
#' @param data Data.frame.
#' @param chr_name Name of column containing chromosome.
#' @param start_name Name of column containing locus start position.
#' @param end_name Name of column containing locus end position.
#' @return GenomicRages object.
#' @importFrom dplyr all_of "%>%"
#' @export

DF2GR <- function(
  data,
  chr_name = "chr",
  start_name = "start",
  end_name = "end"
) {

  if (!methods::is(data, "data.frame")) {
    stop("Input is not a data frame.")
  }
  
  # Rename.
  data <- data %>%
    dplyr::rename(
      chr = all_of({{ chr_name }}),
      start = all_of({{ start_name }}),
      end = all_of({{ end_name }})
    )

  # GRanges object.
  out <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(data$chr),
    ranges = IRanges::IRanges(
      start = data$start,
      end = data$end
    )
  )
  
  # Meta-data columns.
  meta_data_columns <- colnames(data)
  meta_data_columns <- setdiff(
    meta_data_columns,
    c("chr", "start", "end")
  )
  meta_data <- data %>%
    dplyr::select(all_of(meta_data_columns))
  GenomicRanges::mcols(out) <- meta_data
  
  return(out)
}


#' Convert from Genomic Ranges to Data Frame
#'
#' Reverses the action of \code{\link{DF2GR}}.
#'
#' @param gr Genomic Ranges object.
#' @param chr_name Name of column containing chromosome.
#' @param start_name Name of column containing locus start position.
#' @param end_name Name of column containing locus end position.
#' @return Data.frame
#' @importFrom dplyr "%>%"
#' @importFrom methods is
#' @importFrom rlang ":="
#' @export

GR2DF <- function(
  gr,
  chr_name = "chr",
  start_name = "start",
  end_name = "end"
) {

  if (!methods::is(gr, "GRanges")) {
    stop("Input is not a GRanges object.")
  }
  
  width <- NULL
  strand <- NULL
  seqnames <- NULL
  start <- NULL
  end <- NULL
  
  out <- as.data.frame(gr) %>%
    dplyr::select(-width, -strand) %>%
    dplyr::rename(
      !!chr_name := seqnames,
      !!start_name := start,
      !!end_name := end
    )
  return(out)
}
