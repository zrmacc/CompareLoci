#' Overlap class
#'
#' @slot left Input loci on the left.
#' @slot right Input loci on the right.
#' @slot overlaps Results of \code{FindOverlaps(left, right)}.
#' @slot left_overlapping Unique loci on the left overlapping at least 1 locus
#'   on the right.
#' @slot right_overlapped Unique loci on the right overlapped by at least 1
#'   locus on the left.
#' @slot n_left Number of input loci on the left.
#' @slot n_right Number of input loci on the right.
#' @exportClass overlap

setClass(
  Class = "overlap",
  representation = representation(
    left = "data.frame",
    right = "data.frame",
    overlaps = "data.frame",
    left_overlapping = "vector",
    right_overlapped = "vector",
    n_left = "integer",
    n_right = "integer"
  )
)


#' Print Method for Overlap Object.
#'
#' Print method for objects of class \code{overlap}.
#'
#' @param x An object of class \code{overlap}.
#' @param ... Unused.
#' @export

print.overlap <- function(x, ...) {

  cat("Loci on the left overlapping at least 1 locus on the right:\n")
  n_left_overlapping <- length(x@left_overlapping)
  n_left <- x@n_left
  pct <- signif(n_left_overlapping / n_left * 100, digits = 3)
  print(paste0(n_left_overlapping, " of ", n_left, " (", pct, "%)"))
  cat("\n\n")

  cat("Loci on the right overlapped by at least 1 locus on the left:\n")
  n_right_overlapped <- length(x@right_overlapped)
  n_right <- x@n_right
  pct <- signif(n_right_overlapped / n_right * 100, digits = 3)
  print(paste0(n_right_overlapped, " of ", n_right, " (", pct, "%)"))
  cat("\n\n")
}


#' Show Method for Overlap Object
#'
#' @param object An object of class \code{overlap}.
#' @rdname overlap-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "overlap"),
  definition = function(object) {
    print.overlap(x = object)
  }
)
