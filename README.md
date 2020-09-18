# DNAModAnnot

## Introduction
DNAModAnnot is a R package providing a comprehensive toolkit for the genome-wide analysis and annotation of DNA modifications (e.g. 6-methyladenine (6mA)). Its modular architecture allows the analysis of modification detection performed using Pacific Biosciences (PacBio) kineticsTools or Oxford Nanopore Technologies via DeepSignal software.

## Installation
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('Biostrings', 'BSgenome', 'Gviz', 'Logolas'))`

setwd("path/to/package/file/")
install.packages("DNAModAnnot_0.0.0.9015.tar.gz", repos = NULL, type = 'source')
```

You should then be able to load the package into your R session with:

`library(DNAModAnnot)`

## Usage
For detailed instructions, check the package vignette (in 'doc' directory).

## Citation
...
