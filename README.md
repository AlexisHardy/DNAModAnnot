# DNAModAnnot: DNA Modification filtering and Annotation using long-read sequencing data

## Introduction
DNAModAnnot is a R package providing a comprehensive toolkit for the genome-wide analysis and annotation of DNA modifications (e.g. 6-methyladenine (6mA)). Its modular architecture allows the analysis of modification detection performed using Pacific Biosciences (PacBio) kineticsTools or Oxford Nanopore Technologies via DeepSignal software.

## Dependencies

R (>= 4.0.0)

Packages:
Biostrings (>= 2.10.0), GenomicRanges (>= 1.38.0), BSgenome (>= 1.28.0), Biobase (> 2.1.0), 
BiocGenerics (>= 0.34.0), GenomeInfoDb (>= 1.14.0), Gviz (>= 1.29.1), IRanges (>= 2.20.0), 
Logolas (>= 1.3.1), S4Vectors (>= 0.24.0), data.table (>= 1.13.0), rtracklayer (>= 1.30.0), 
Rsamtools (>= 2.0.0)

## Installation

First, install required packages using BiocManager:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('Biostrings', 'BSgenome', 'Gviz', 'Logolas'))
```

Then, you can install DNAModAnnot using the tar.gz file from GitHub repository:
```
setwd("path/to/package/file/")
install.packages("DNAModAnnot_0.0.0.9015.tar.gz", repos = NULL, type = 'source')
```

Or you can directly install from GitHub using devtools package:
```
install.packages("devtools")
library(devtools)
install_github("AlexisHardy/DNAModAnnot")
```

You should then be able to load the package into your R session with:
`library(DNAModAnnot)`

## Usage
For detailed instructions, check the package vignette (in 'doc' directory).

## Citation
...
