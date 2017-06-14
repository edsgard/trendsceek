#' Single-cell RNA-seq data from the Scialdone et al paper
#'
#' A dataset containing read counts for the subset of 481 cells annotated as "cluster 3" in the paper by Scialdone et al.
#' 
#' @format A list with two elements, "counts" and "meta". Counts is a data-frame containing read counts for 41,388 genes and 481 cells. Meta is a data-frame containing annotations about all 481 cells.
#' 
#' @source \itemize{
#' \item RNA-seq read counts \url{http://gastrulation.stemcells.cam.ac.uk/data/counts.gz}
#' \item{RNA-seq meta-info}{\url{http://gastrulation.stemcells.cam.ac.uk/data/metadata.txt}}
#' }
#' @usage data(scialdone)
"scialdone"
