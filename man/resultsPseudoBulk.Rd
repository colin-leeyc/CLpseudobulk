% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk_functions.R
\name{resultsPseudoBulk}
\alias{resultsPseudoBulk}
\title{This is an internal function used by iteratePseudoBulk}
\usage{
resultsPseudoBulk(
  sce_int,
  cellIDs_int = NULL,
  obs_int,
  contrasts_int,
  replicates_int = NULL,
  alpha_int = 0.05,
  genes_int = NULL
)
}
\value{
A list of LogFoldChanges, adjusted P values, and p values for genes of interest
}
\description{
This is an internal function used by iteratePseudoBulk
}
