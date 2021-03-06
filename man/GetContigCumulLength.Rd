% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeQual.R
\name{GetContigCumulLength}
\alias{GetContigCumulLength}
\title{GetContigCumulLength Function (SeQual)}
\usage{
GetContigCumulLength(dnastringsetGenome)
}
\arguments{
\item{dnastringsetGenome}{A DNAStringSet object containing the sequence for each contig.}
}
\value{
A dataframe:
\itemize{
  \item contig_names: The names of each contig.
  \item Mbp_length: The size of each contig (in Megabase pairs (Mbp)).
  \item cumsum_Mbp_length: The cumulative size (from largest to smallest contig) for each contig (in Megabase pairs (Mbp)).
}
}
\description{
Return a dataframe with the length (and cumulative length) of contigs (ordered from largest to smallest contig) provided in genome assembly.
}
\examples{
myGenome <- Biostrings::readDNAStringSet(system.file(
  package = "DNAModAnnot", "extdata",
  "ptetraurelia_mac_51_sca171819.fa"
))

myContig_cumul_len_t <- GetContigCumulLength(dnastringsetGenome = myGenome)

myContig_cumul_len_t
}
\keyword{GetContigCumulLength}
