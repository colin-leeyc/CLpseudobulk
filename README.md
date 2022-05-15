# CLpseudobulk
Pseudobulk analysis for scRNAseq by creating artifical replicates by random sampling.
Subsequently able to iterate random sampling (bootstrapping) and workflow to robustly test the differential expression between conditions.

Please install the following packages in R prior to use:

SingleCellExperiment,
DESeq2,
forcats,
dplyr,
tibble.

To install: 
devtools::install_github("colin-leeyc/CLpseudobulk")

