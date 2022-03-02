# Compare Loci

Zachary McCaw <br>
Updated: 2022-03-02



### Description

Simple application of the [GRanges](https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html) package to compare two sets of loci.

## Installation


```r
devtools::install_github(repo = "zrmacc/CompareLoci")
```

## Example


```r
library(CompareLoci)

# Sets of loci to compare. 
left <- data.frame(
  chr = c(1, 1, 1),
  start = c(1, 200, 400),
  end = c(100, 300, 500)
)

right <- data.frame(
  chr = c(1, 1, 1, 1),
  start = c(20, 50, 150, 440),
  end = c(25, 75, 180, 450)
)

# Overlap statistics.
overlap_stats <- CalcOverlapStats(left, right)
show(overlap_stats)
```

```
## Loci on the left overlapping at least 1 locus on the right:
## [1] "2 of 3 (66.7%)"
## 
## 
## Loci on the right overlapped by at least 1 locus on the left:
## [1] "3 of 4 (75%)"
```

```r
# Loci on the left overlapping at least 1 locus on the right:
show(overlap_stats@left_overlapping)
```

```
## [1] 1 3
```

```r
# Loci on the right overlapping at least 1:
show(overlap_stats@right_overlapped)
```

```
## [1] 1 2 4
```
