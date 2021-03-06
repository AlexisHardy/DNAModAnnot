% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeQual.R
\name{GetSeqPctByContig}
\alias{GetSeqPctByContig}
\title{GetSeqPctByContig Function (SeQual)}
\usage{
GetSeqPctByContig(gposPacBioCSV, grangesGenome)
}
\arguments{
\item{gposPacBioCSV}{An UnStitched GPos object containing PacBio CSV data to be analysed.}

\item{grangesGenome}{A GRanges object containing the width of each contig.}
}
\value{
A list composed of 3 dataframes: 1 dataframe by strand and 1 dataframe with both strands. In each dataframe:
\itemize{
  \item refName: The names of each contig.
  \item strand: The strand of each contig.
  \item width: The width of each contig.
  \item nb_sequenced: The number of bases sequenced by strand for each contig.
  \item seqPct: The percentage of bases sequenced for each strand for each contig (percentage of sequencing).
}
}
\description{
Return a list with the percentage of sequencing by strand for all scaffolds of genome assembly provided.
This function is not adapted for data from DeepSignal.
}
\examples{
myGenome <- Biostrings::readDNAStringSet(system.file(
  package = "DNAModAnnot", "extdata",
  "ptetraurelia_mac_51_sca171819.fa"
))
myGrangesGenome <- GetGenomeGRanges(myGenome)

# Preparing a gposPacBioCSV dataset
myGposPacBioCSV <-
  ImportPacBioCSV(
    cPacBioCSVPath = system.file(
      package = "DNAModAnnot", "extdata",
      "ptetraurelia.bases.sca171819.csv"
    ),
    cSelectColumnsToExtract = c(
      "refName", "tpl", "strand", "base",
      "score", "ipdRatio", "coverage"
    ),
    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
    cContigToBeAnalyzed = names(myGenome)
  )

myPct_seq_csv <- GetSeqPctByContig(myGposPacBioCSV, grangesGenome = myGrangesGenome)
myPct_seq_csv
}
\keyword{GetSeqPctByContig}
